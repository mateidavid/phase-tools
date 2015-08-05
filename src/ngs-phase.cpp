#include <iostream>
#include <sstream>
#include <cstdlib>
#include <time.h>
#include <string>
#include <set>
#include <vector>

#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <tclap/CmdLine.h>
#include "version.h"

#include "zstr.hpp"
#include "logger.hpp"
#include "RC_Sequence.hpp"
#include "Het_Variation.hpp"
#include "Phased_Set.hpp"
#include "Phased_Set_Connection.hpp"
#include "Cigar.hpp"
#include "Mapping.hpp"
#include "overlapper.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "missing_version"
#endif

using namespace std;

typedef rc_sequence::Sequence< string > DNA_Sequence;

namespace global
{
    using namespace TCLAP;

    string description =
        "Given a VCF file of variants and one or more BAM files with NGS reads, "
        "compute phased genotypes directly supported by the reads.";
    CmdLine cmd_parser(description, ' ', PACKAGE_VERSION);
    //
    // general parameters
    //
    MultiArg< string > log_level("d", "log-level", "Log level.", false, "string", cmd_parser);
    ValueArg< long > seed("", "seed", "Random seed (-1: use time).", false, -1, "int", cmd_parser);
    //
    // io-related parameters
    //
    ValueArg< string > ref_fn("", "ref", "Reference (Fasta) file.", true, "", "file", cmd_parser);
    ValueArg< string > var_fn("", "var", "Variants (VCF) file.", true, "", "file", cmd_parser);
    MultiArg< string > map_fn("", "map", "Mappings (BAM) file.", true, "file", cmd_parser);
    ValueArg< string > out_var_fn("", "out", "Output variations file", true, "", "file", cmd_parser);
    ValueArg< string > sample_id("", "sample", "Sample Id.", true, "", "string", cmd_parser);
    //
    // other parameters
    //
    ValueArg< string > chr("", "chr", "Process single chromosome.", false, "", "chr", cmd_parser);
    MultiArg< string > skip_chr("", "skip-chr", "Skip chromosome.", false, "chr", cmd_parser);
    ValueArg< int > flank_len("", "flank_len", "Flank length [20].", false, 20, "int", cmd_parser);
    ValueArg< double > max_discordance("", "max-discordance", "Maximum discordance for connecting phase sets [.2].", false, .2, "float", cmd_parser);
    ValueArg< string > gt_tag("", "gt-tag", "GT tag.", false, "GT_ngs", "string", cmd_parser);
    ValueArg< string > ps_tag("", "ps-tag", "PS tag.", false, "PS_ngs", "string", cmd_parser);

    // reference
    faidx_t * faidx_p;
    // phased chromosomes
    set< string > skip_chr_s;
    // variations
    map< string, deque< Het_Variation > > var_m;
    map< string, set< const Het_Variation *, Het_Variation_Ptr_Comp > > het_m;
    // list of chromosomes in order of input file
    list< string > chrom_l;
    // list of genotypes observed by a fragment
    map< string, map< const Het_Variation *, bool > > frag_store_m;
    // output structures
    vcfFile * of_p;
    bcf_hdr_t * ohdr_p;

    // counts

    size_t num_in_map_total;
    size_t num_in_map_nonprimary;
    size_t num_in_map_unpaired;
    size_t num_in_map_unpaired_unmapped;
    size_t num_in_map_unpaired_mapped;
    size_t num_in_map_paired;
    size_t num_in_map_paired_unmapped;
    size_t num_in_map_paired_mapped;
    size_t num_in_map_paired_both_unmapped;
    size_t num_in_map_paired_both_mapped;
    size_t num_in_map_paired_both_mapped_diff_chr;

    size_t num_out_map_unpaired_total;
    size_t num_out_map_unpaired_inconclusive;
    size_t num_out_map_unpaired_random;
    size_t num_out_map_unpaired_phased;

    size_t num_out_frag_total;
    size_t num_out_frag_conflicting;
    size_t num_out_frag_inconclusive;
    size_t num_out_frag_random;
    size_t num_out_frag_single_sided;
    size_t num_out_frag_concordant;
    size_t num_out_frag_missing_unmapped;
    size_t num_out_frag_missing_mapped;

    size_t num_out_map_total;
    size_t num_out_map_skip_chr;
    size_t indel_phasing_neither_flank;

} // namespace global

void set_skip_chr_s()
{
    for (const auto & s : global::skip_chr.get())
    {
        LOG("main", info) << "skipping chromosome [" << s << "]" << endl;
        global::skip_chr_s.insert(s);
    }
}

void load_reference_index()
{
    LOG("main", info) << "load_reference_index(): start" << endl;
    global::faidx_p = fai_load(global::ref_fn.get().c_str());
    LOG("main", info) << "load_reference_index(): end" << endl;
}

void load_variations()
{
    LOG("main", info) << "load_variations(): start" << endl;
    int ret;
    // open vcf file, initialize structs
    auto f_p = bcf_open(global::var_fn.get().c_str(), "r");
    auto hdr_p = bcf_hdr_read(f_p);
    global::ohdr_p = bcf_hdr_init("w");
    auto rec_p = bcf_init1();
    int * dat = nullptr;
    int dat_size = 0;
    // restrict processing to given sample_id
    ret = bcf_hdr_set_samples(hdr_p, global::sample_id.get().c_str(), 0);
    if (ret != 0)
    {
        cerr << "bcf_hdr_set_samples() error: status=" << ret << endl;
        exit(1);
    }

    // set up output header (copy contigs)
    bcf_hdr_append(global::ohdr_p, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(global::ohdr_p, "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase Set\">");
    int n_seqs;
    auto seq_names = bcf_hdr_seqnames(hdr_p, &n_seqs);
    for (int i = 0; i < n_seqs; ++i)
    {
        auto hrec_p = bcf_hdr_get_hrec(hdr_p, BCF_HL_CTG, "ID", seq_names[i], nullptr);
        bcf_hdr_add_hrec(global::ohdr_p, bcf_hrec_dup(hrec_p));
    }
    free(seq_names);
    bcf_hdr_add_sample(global::ohdr_p, global::sample_id.get().c_str());
    bcf_hdr_add_sample(global::ohdr_p, nullptr); //?!?!

    // main loop
    while (bcf_read1(f_p, hdr_p, rec_p) >= 0)
    {
        bcf_unpack(rec_p, BCF_UN_ALL);
        Het_Variation v(hdr_p, rec_p, dat, dat_size, false);
        string chrom = v.chr_name();
        if (global::skip_chr_s.count(chrom) > 0)
        {
            continue;
        }
        if (global::chrom_l.empty() or chrom != global::chrom_l.back())
        {
            global::chrom_l.push_back(chrom);
        }
        if (v.is_het())
        {
            v.reset_phase();
            v.load_flanks(global::faidx_p, global::flank_len);
        }
        LOG("variations", debug) << v << endl;
        global::var_m[chrom].push_back(move(v));
        const Het_Variation & v_ref = global::var_m[chrom].back();
        if (v_ref.is_het())
        {
            auto res = global::het_m[chrom].insert(&v_ref);
            if (not res.second)
            {
                cerr << "duplicate variation:\n" << v_ref << endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    // destroy structures and close file
    if (dat) free(dat);
    //delete [] dat;
    bcf_destroy1(rec_p);
    bcf_hdr_destroy(hdr_p);
    ret = bcf_close(f_p);
    if (ret != 0)
    {
        cerr << "bcf_close() error: status=" << ret << endl;
        exit(ret);
    }
    LOG("main", info) << "load_variations(): end" << endl;
}

void output_variations_chromosome(const string & chrom)
{
    int odat[2];
    for (const auto & v : global::var_m[chrom])
    {
        assert(v.chr_name() == chrom);
        auto orec_p = bcf_init1();
        // set chr and pos
        orec_p->rid = bcf_hdr_name2id(global::ohdr_p, chrom.c_str());
        orec_p->pos = v.rf_start();
        // set alleles
        ostringstream oss;
        for (size_t i = 0; i < v.n_alleles(); ++i)
        {
            oss << (i > 0? "," : "") << v.allele_seq(i);
        }
        bcf_update_alleles_str(global::ohdr_p, orec_p, oss.str().c_str());
        if (v.is_het())
        {
            // set GT
            odat[0] = bcf_gt_phased(v.gt(v.ps_phase));
            odat[1] = bcf_gt_phased(v.gt(1 - v.ps_phase));
            bcf_update_genotypes(global::ohdr_p, orec_p, odat, 2);
            // set PS
            bcf_update_format_int32(global::ohdr_p, orec_p, "PS", &v.ps_start_1, 1);
        }
        else
        {
            // set GT
            bcf_update_genotypes(global::ohdr_p, orec_p, v.raw_gt().data(), v.raw_gt().size());
        }
        bcf_write1(global::of_p, global::ohdr_p, orec_p);
        bcf_destroy1(orec_p);
    }

    /*
    // transform output file
    while (bcf_read1(if_p, ihdr_p, irec_p) >= 0)
    {
        bcf_unpack(irec_p, BCF_UN_ALL);
        Het_Variation v2(ihdr_p, irec_p, dat, dat_size, false);
        if (not v2.is_valid()
            or (not global::chr.get().empty() and v2.chr_name() != global::chr.get()))
        {
            continue;
        }
        auto it = global::var_m[v2.chr_name()].find(Het_Variation(v2.rf_start()));
        if (it == global::var_m[v2.chr_name()].end())
        {
            LOG("main", warning) << "did not find the following variation loaded during output:" << endl
                                 << v2 << endl;
            continue;
        }
        const Het_Variation & v = *it;
        // update output record
        auto orec_p = bcf_init1();
        orec_p->rid = bcf_hdr_name2id(ohdr_p, v.chr_name().c_str());
        orec_p->pos = irec_p->pos;
        ostringstream oss;
        for (size_t i = 0; i < v.n_alleles(); ++i)
        {
            oss << (i > 0? "," : "") << v.allele_seq(i);
        }
        bcf_update_alleles_str(ohdr_p, orec_p, oss.str().c_str());
        // update GT
        odat[0] = bcf_gt_phased(v.gt(v.ps_phase));
        odat[1] = bcf_gt_phased(v.gt(1 - v.ps_phase));
        bcf_update_genotypes(ohdr_p, orec_p, odat, 2);
        // update PS
        bcf_update_format_int32(ohdr_p, orec_p, "PS", &v.ps_start_1, 1);
        bcf_write1(of_p, ohdr_p, orec_p);
        bcf_destroy1(orec_p);
    }
    */
}

void output_variations()
{
    static const vector< pair< string, string > > extensions{ {".vcf.gz", "wz"}, {".bcf", "wb"}, {".vcf", "w"} };
    string output_mode = "w";
    for (const auto & p : extensions)
    {
        if (global::out_var_fn.get().size() >= p.first.size() + 1
            and global::out_var_fn.get().substr(global::out_var_fn.get().size() - p.first.size()) == p.first)
        {
            output_mode = p.second;
            break;
        }
    }
    global::of_p = bcf_open(global::out_var_fn.get().c_str(), output_mode.c_str());
    // write out header
    bcf_hdr_write(global::of_p, global::ohdr_p);
    // process chromosomes
    if (global::chr.get().empty())
    {
        for (const auto & chrom : global::chrom_l)
        {
            output_variations_chromosome(chrom);
        }
    }
    else
    {
        output_variations_chromosome(global::chr);
    }
    // destroy structures and close file
    bcf_hdr_destroy(global::ohdr_p);
    bcf_close(global::of_p);
}

void print_variations(ostream& os)
{
    for (const auto & p : global::var_m)
    {
        for (const auto & v : p.second)
        {
            os << v;
        }
    }
}

int get_phase_snp(const Mapping & m, const Het_Variation & v)
{
    assert(v.is_snp());
    assert(m.rf_start() <= v.rf_start() and v.rf_end() <= m.rf_start() + m.rf_len());
    // compute the mapped range, on query, of the single reference position
    auto rg = m.cigar().mapped_range(make_pair(v.rf_start() - m.rf_start(), v.rf_start() - m.rf_start() + 1), true, true);
    if (rg.second != rg.first + 1)
    {
        // the changed reference base is mapped to something other than a match; we give up
        return -1;
    }
    else
    {
        // the changed reference base is mapped to seq[rg.first]
        assert(0 <= rg.first and rg.second <= static_cast< int >(m.seq().size()));
        if (m.seq()[rg.first] == v.allele_seq(v.gt(0))[0])
            return 0;
        else if (m.seq()[rg.first] == v.allele_seq(v.gt(1))[0])
            return 1;
        else
            return -1;
    }
}

int get_phase_indel(const Mapping & m, const Het_Variation & v)
{
    // We first estimate the mapping of the query to the flanking regions of the
    // indel.
    // delta = estimate of the max indel size in flank mapping
    // We use delta to compute a band limit for the DP.
    pair< int, int > rf_rg[2];
    pair< int, int > qr_rg[2];
    int delta[2];
    assert(v.rf_start() <= m.rf_end());
    assert(m.rf_start() <= v.rf_end());
    rf_rg[0] = make_pair(max(v.rf_start() - global::flank_len - m.rf_start(), 0),
                         max(v.rf_start() - m.rf_start(), 0));
    rf_rg[1] = make_pair(min(v.rf_end() - m.rf_start(), m.rf_len()),
                         min(v.rf_end() + global::flank_len - m.rf_start(), m.rf_len()));
    qr_rg[0] = m.cigar().mapped_range(rf_rg[0], true, false);
    qr_rg[1] = m.cigar().mapped_range(rf_rg[1], true, false);
    delta[0] = abs((qr_rg[0].second - qr_rg[0].first) - (rf_rg[0].second - rf_rg[0].first));
    delta[1] = abs((qr_rg[1].second - qr_rg[1].first) - (rf_rg[1].second - rf_rg[1].first));

    DNA_Sequence rf_allele[2];
    rf_allele[0] = v.flank_seq(0) + v.allele_seq(v.gt(0)) + v.flank_seq(1);
    rf_allele[1] = v.flank_seq(0) + v.allele_seq(v.gt(1)) + v.flank_seq(1);
    int max_allele_seq_size = max(v.allele_seq(v.gt(0)).size(), v.allele_seq(v.gt(1)).size());
    DNA_Sequence qr;
    bool anchor_head = true;
    int band_width;
    if (rf_rg[0].second - rf_rg[0].first == global::flank_len
        and (rf_rg[1].second - rf_rg[1].first < global::flank_len
             or delta[0] < delta[1]))
    {
        // use 5p flank
        int qr_end = qr_rg[0].second;
        qr_end = min(qr_end + max_allele_seq_size + global::flank_len, static_cast< int >(m.seq().size()));
        qr = m.seq().substr(qr_rg[0].first, qr_end - qr_rg[0].first);
        anchor_head = true;
        band_width = delta[0];
    }
    else if (rf_rg[1].second - rf_rg[1].first == global::flank_len)
    {
        // use 3p flank
        int qr_start = qr_rg[1].first;
        qr_start = max(qr_start - max_allele_seq_size - global::flank_len, 0);
        qr = m.seq().substr(qr_start, qr_rg[1].second - qr_start);
        anchor_head = false;
        band_width = delta[1];
    }
    else
    {
        // neither flank is fully captured by the mapping, we give up
        ++global::indel_phasing_neither_flank;
        return -1;
    }
    SequenceOverlap res;
    int score[2];
    res = Overlapper::extendMatch(rf_allele[0], qr,
                                  (anchor_head? 0 : rf_allele[0].size() - 1),
                                  (anchor_head? 0 : qr.size() - 1),
                                  band_width);
    score[0] = res.score;
    res = Overlapper::extendMatch(rf_allele[1], qr,
                                  (anchor_head? 0 : rf_allele[1].size() - 1),
                                  (anchor_head? 0 : qr.size() - 1),
                                  band_width);
    score[1] = res.score;
    return (score[0] == score[1]? -1 : (score[0] > score[1]? 0 : 1));
}

int get_phase(const Mapping & m, const Het_Variation & v)
{
    // for now, only deal with SNPs
    int res = v.is_snp()
        ? get_phase_snp(m, v)
        : get_phase_indel(m, v);
    // update counts
    ++v.frag_total;
    if (res >= 0) ++v.frag_supp_allele[res];
    return res;
}

/**
 * Process one mapping.
 */
void process_mapping(const Mapping & m)
{
    if (global::skip_chr_s.count(m.chr_name()))
    {
        ++global::num_out_map_skip_chr;
        return;
    }
    if (not global::het_m.count(m.chr_name()))
    {
        // no hets on this chr
        return;
    }
    if (not m.is_mapped())
    {
        return;
    }
    // compute range of variations spanned by this mapping
    Het_Variation het_start(m.rf_start());
    Het_Variation het_end(m.rf_end());
    auto it_start = global::het_m.at(m.chr_name()).lower_bound(&het_start);
    while (it_start != global::het_m.at(m.chr_name()).begin()
           and (*prev(it_start))->rf_end() > m.rf_start())
    {
        --it_start;
    }
    auto it_end = global::het_m.at(m.chr_name()).lower_bound(&het_end);
    if (it_start == it_end)
    {
        // not spanning any hets
        return;
    }
    // mapping spans at least 1 het
    auto & fm = global::frag_store_m[m.query_name()];
    for (auto it = it_start; it != it_end; ++it)
    {
        const Het_Variation & v = **it;
        int phase = get_phase(m, v);
        if (phase >= 0)
        {
            LOG("phasing", debug1)
                << "fragment [" << m.query_name()
                << "] observes variation [" << v.chr_name() << ":" << v.rf_start() + 1 << "] phase ["
                << phase << "]" << endl;
            if (fm.count(&v))
            {
                if (phase != fm[&v])
                {
                    LOG("main", warning)
                        << "fragment [" << m.query_name() << "] assigned to both phases of variation ["
                        << v.chr_name() << ":" << v.rf_start() << "]" << endl;
                    fm.erase(&v);
                }
            }
            else
            {
                fm[&v] = phase;
            }
        }
    }
}

void process_mappings(const string & fn, const string & chr)
{
    htsFile * map_file_p = hts_open(fn.c_str(), "r");
    bam_hdr_t * map_hdr_p = sam_hdr_read(map_file_p);
    bam1_t * rec_p = bam_init1();
    std::function< int(void) > get_next_record;

    // traverse given chromosome only
    hts_idx_t * map_idx_p = sam_index_load(map_file_p, fn.c_str());
    if (not map_idx_p)
    {
        cerr << "could not load BAM index for [" << fn << "]" << endl;
        exit(EXIT_FAILURE);
    }
    int tid = bam_name2id(map_hdr_p, chr.c_str());
    if (tid < 0)
    {
        LOG("main", warning) << "chromosome [" << chr << "] not found in BAM header for ["
                             << fn << "]" << endl;
    }
    hts_itr_t * map_itr_p = sam_itr_queryi(map_idx_p, tid, 0, numeric_limits< int >::max());
    get_next_record = [&] () { return sam_itr_next(map_file_p, map_itr_p, rec_p); };

    while (get_next_record() >= 0)
    {
        Mapping m(map_hdr_p, rec_p);
        process_mapping(m);
    }

    hts_itr_destroy(map_itr_p);
    hts_idx_destroy(map_idx_p);
    bam_destroy1(rec_p);
    bam_hdr_destroy(map_hdr_p);
    hts_close(map_file_p);
}

void process_chromosome(const string & chr)
{
    assert(global::frag_store_m.empty());
    if (global::skip_chr_s.count(chr)) return;

    // load mappings evidence from all bams
    for (const auto & fn : global::map_fn)
    {
        process_mappings(fn, chr);
    }

    // create initial phased sets as single hets
    vector< Phased_Set > phased_set_v(global::het_m.at(chr).size());
    size_t i = 0;
    for (const auto & v_p : global::het_m.at(chr))
    {
        phased_set_v[i].het_set().insert(make_pair(v_p, false));
        v_p->phased_set_ptr = &phased_set_v[i];
        ++i;
    }
    // initialize collection of phased set connections
    Connection_Store conn_store;
    Connection_Map conn_map;
    for (const auto & p : global::frag_store_m)
    {
        const auto & fm = p.second;
        for (auto it1 = fm.begin(); it1 != fm.end(); ++it1)
            for (auto it2 = next(it1); it2 != fm.end(); ++it2)
            {
                assert(it1->first != it2->first);
                const Het_Variation & v1 = *it1->first;
                bool v1_phase = it1->second;
                const Het_Variation & v2 = *it2->first;
                bool v2_phase = it2->second;
                // fragment connects v1:v1_phase to v2:v2_phase
                Phased_Set & ps1 = *v1.phased_set_ptr;
                bool ps1_phase = ps1.het_set().find(make_pair(&v1, false)) == ps1.het_set().end();
                assert(ps1.het_set().find(make_pair(&v1, ps1_phase)) != ps1.het_set().end());
                Phased_Set & ps2 = *v2.phased_set_ptr;
                bool ps2_phase = ps2.het_set().find(make_pair(&v2, false)) == ps2.het_set().end();
                assert(ps2.het_set().find(make_pair(&v2, ps2_phase)) != ps2.het_set().end());
                // thus ps1:(v1_phase + ps1_phase) to ps2:(v2_phase + ps2_phase)
                ps1_phase = (v1_phase + ps1_phase) % 2;
                ps2_phase = (v2_phase + ps2_phase) % 2;
                // now, ps1:ps1_phase to ps2:ps2_phase
                LOG("phasing", debug1)
                    << "fragment [" << p.first << "] connects ["
                    << ps1 << ":" << ps1_phase << "] with ["
                    << ps2 << ":" << ps2_phase << "]" << endl;
                Phased_Set_Connection * cp = nullptr;
                if (conn_map[&ps1].count(&ps2) == 0)
                {
                    conn_store.emplace_back(&ps1, &ps2);
                    cp = &conn_store.back();
                    conn_map[&ps1][&ps2] = cp;
                    conn_map[&ps2][&ps1] = cp;
                }
                else
                {
                    cp = conn_map[&ps1][&ps2];
                    assert(cp == conn_map[&ps2][&ps1]);
                }
                cp->increment_count(ps1_phase, ps2_phase);
            }
    }
    // build connection set, sorted by discordance
    Connection_Set conn_set;
    for (auto & c : conn_store)
    {
        conn_set.insert(c);
    }

    // main greedy loop
    while (not conn_set.empty())
    {
        // pick connection with smallest discordance
        Phased_Set_Connection & c = *conn_set.begin();
        LOG("phasing", debug1) << "connection with smallest discordance: " << c << endl;
        // stop if smallest discordance is too large
        if (c.discordance() > global::max_discordance) break;
        // connect phased sets according to the current connection
        Phased_Set * ps1_ptr = const_cast< Phased_Set * >(c.ps_ptr(0));
        Phased_Set * ps2_ptr = const_cast< Phased_Set * >(c.ps_ptr(1));
        if (ps1_ptr->het_set().size() < ps2_ptr->het_set().size()) swap(ps1_ptr, ps2_ptr);
        bool rel_phase = c.rel_phase();
        // merge ps2 into ps1
        for (const auto & p : ps2_ptr->het_set())
        {
            ps1_ptr->het_set().insert(make_pair(p.first, p.second != rel_phase));
            p.first->phased_set_ptr = ps1_ptr;
        }
        ps2_ptr->het_set().clear();
        // remove connection ps1 to ps2
        conn_set.erase(conn_set.iterator_to(c));
        conn_map[ps1_ptr].erase(ps2_ptr);
        conn_map[ps2_ptr].erase(ps1_ptr);
        // transform connections to ps2 into connections to ps1
        for (const auto & p : conn_map[ps2_ptr])
        {
            Phased_Set * ps3_ptr = const_cast< Phased_Set * >(p.first);
            Phased_Set_Connection * old_cp = p.second;
            conn_set.erase(conn_set.iterator_to(*old_cp));
            Phased_Set_Connection * new_cp = nullptr;
            if (conn_map[ps1_ptr].count(ps3_ptr) == 0)
            {
                // simply move connection
                old_cp->replace_ps_ptr(ps2_ptr, ps1_ptr, rel_phase);
                conn_set.insert(*old_cp);
                new_cp = old_cp;
            }
            else
            {
                // merge connections
                Phased_Set_Connection * direct_cp = conn_map[ps1_ptr][ps3_ptr];
                conn_set.erase(conn_set.iterator_to(*direct_cp));
                direct_cp->merge(*old_cp, ps3_ptr, rel_phase);
                conn_set.insert(*direct_cp);
                new_cp = direct_cp;
            }
            assert(conn_map[ps3_ptr].count(ps2_ptr) > 0);
            conn_map[ps3_ptr].erase(ps2_ptr);
            conn_map[ps1_ptr][ps3_ptr] = new_cp;
            conn_map[ps3_ptr][ps1_ptr] = new_cp;
        }
        conn_map.erase(ps2_ptr);
        LOG("phasing", debug1) << "merged phased set: " << *ps1_ptr << endl;
    }
    if (Logger::get_facility_level("phasing") >= Logger::level::debug1)
    {
        for (const auto & c : conn_set)
        {
            LOG("phasing", debug1) << "discarded connection: " << c << endl;
        }
    }

    // output
    for (const auto & v_p : global::het_m[chr])
    {
        Phased_Set & ps = *v_p->phased_set_ptr;
        v_p->ps_phase = ps.het_set().find(make_pair(v_p, false)) == ps.het_set().end();
        v_p->ps_start_1 = ps.het_set().begin()->first->rf_start() + 1;
    }

    // clear intrusive structure (the rest are cleared automatically)
    conn_set.clear();

    global::frag_store_m.clear();
}

void print_stats(ostream & os)
{
    os
        << "input_mappings: " << global::num_in_map_total << endl
        << "..nonprimary: " << global::num_in_map_nonprimary << endl
        << "..unpaired: " << global::num_in_map_unpaired << endl
        << "....unmapped: " << global::num_in_map_unpaired_unmapped << endl
        << "....mapped: " << global::num_in_map_unpaired_mapped << endl
        << "..paired: " << global::num_in_map_paired << endl
        << "....unmapped: " << global::num_in_map_paired_unmapped << endl
        << "....mapped: " << global::num_in_map_paired_mapped << endl
        << "....both_unmapped: " << global::num_in_map_paired_both_unmapped << endl
        << "....both_mapped: " << global::num_in_map_paired_both_mapped << endl
        << "......diff_chr: " << global::num_in_map_paired_both_mapped_diff_chr << endl

        << "output_mappings: " << global::num_out_map_total << endl
        << "..unpaired: " << global::num_out_map_unpaired_total << endl
        << "....inconclusive: " << global::num_out_map_unpaired_inconclusive << endl
        << "....random: " << global::num_out_map_unpaired_random << endl
        << "....phased: " << global::num_out_map_unpaired_phased << endl

        << "..paired_fragments: " << global::num_out_frag_total << endl
        << "....conflicting: " << global::num_out_frag_conflicting << endl
        << "....inconclusive: " << global::num_out_frag_inconclusive << endl
        << "....random: " << global::num_out_frag_random << endl
        << "....single_sided: " << global::num_out_frag_single_sided << endl
        << "....concordant: " << global::num_out_frag_concordant << endl
        << "....missing_unmapped: " << global::num_out_frag_missing_unmapped << endl
        << "....missing_mapped: " << global::num_out_frag_missing_mapped << endl

        << "..skip_chr: " << global::num_out_map_skip_chr << endl
        ;
}

void print_var_stats(ostream & os)
{
    if (global::chr.get().empty())
    {
        for (const auto & s : global::het_m)
        {
            for(const auto & v_p : s.second)
            {
                os << *v_p << endl;
            }
        }
    }
    else
    {
        if (global::het_m.count(global::chr))
        {
            for (const auto & v_p : global::het_m.at(global::chr))
            {
                os << *v_p << endl;
            }
        }
    }
}

void real_main()
{
    set_skip_chr_s();
    load_reference_index();
    load_variations();
    //print_variations(clog);
    if (global::chr.get().empty())
    {
        for (const auto & chrom : global::chrom_l)
        {
            process_chromosome(chrom);
        }
    }
    else
    {
        process_chromosome(global::chr);
    }
    output_variations();
    // cleanup
    fai_destroy(global::faidx_p);
}

int main(int argc, char * argv[])
{
    global::cmd_parser.parse(argc, argv);
    // set log levels
    Logger::set_levels_from_options(global::log_level, &clog);
    // print options
    LOG("main", info) << "program: " << global::cmd_parser.getProgramName() << endl;
    LOG("main", info) << "version: " << global::cmd_parser.getVersion() << endl;
    LOG("main", info) << "args: " << global::cmd_parser.getOrigArgv() << endl;
    // set random seed
    if (global::seed < 0)
    {
        global::seed.get() = time(nullptr);
        LOG("main", info) << "seed: " << global::seed << endl;
    }
    srand48(global::seed);
    // real main
    real_main();
}
