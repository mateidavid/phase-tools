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

#include "zstr.hpp"
#include "logger.hpp"
#include "RC_Sequence.hpp"
#include "Het_Variation.hpp"
#include "Cigar.hpp"
#include "Mapping.hpp"
#include "Phaser.hpp"
#include "version.hpp"

using namespace std;

namespace global
{
    using namespace TCLAP;

    string description =
        "Given a BAM file and a partly phased VCF file, split the mappings in the BAM file "
        "into several other BAM files. Specifically, for every phase set (PS) defined in the VCF file, "
        "and for every autosome (1-22), there will be 2 BAM output files (1 per haplotype). "
        "In addition, there will be 2 more overflow BAM output files for reads not mapped to autosomes.";
    CmdLine cmd_parser(description, ' ', package_version);
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
    ValueArg< string > map_fn("", "map", "Mappings (BAM) file.", true, "", "file", cmd_parser);
    ValueArg< string > out_map_fn("", "out", "Output file prefix. Files names will be <prefix>.<chr>.[01].bam", true, "", "file", cmd_parser);
    ValueArg< string > sample_id("", "sample", "Sample Id.", true, "", "string", cmd_parser);
    ValueArg< string > var_stats_fn("", "var-stats", "Variations stats file.", false, "", "file", cmd_parser);
    //
    // other parameters
    //
    ValueArg< string > chr("", "chr", "Process single chromosome.", false, "", "chr", cmd_parser);
    MultiArg< string > chr_phase("", "chr-phase", "Assign all mappings to this chr to given phase.", false, "chr:[01]", cmd_parser);
    ValueArg< int > flank_len("", "flank_len", "Flank length [20].", false, 20, "int", cmd_parser);
    ValueArg< string > decision_fn("", "decision", "File to log decision.", false, "", "file", cmd_parser);

    // reference
    faidx_t * faidx_p;
    // phased chromosomes
    map< string, int > chr_phase_m;
    // variations
    map< string, set< Het_Variation > > var_m;
    // input mappings
    htsFile * map_file_p;
    bam_hdr_t * map_hdr_p;
    hts_idx_t * map_idx_p;
    // output mappings
    htsFile * crt_out_map_file_p[2];
    // mate pair decision
    map< string, tuple< const Mapping *, int, vector< pair< const Het_Variation *, int > > > > mp_store_m;
    // decision log
    strict_fstream::ofstream decision_ofs;

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
    size_t num_out_frag_hets_neither_side;
    size_t num_out_frag_hets_single_side;
    size_t num_out_frag_hets_single_side_inconclusive;
    size_t num_out_frag_hets_both_sides;
    size_t num_out_frag_hets_both_sides_inconclusive;
    size_t num_out_frag_random;
    size_t num_out_frag_missing_mapped;
    size_t num_out_frag_missing_unmapped;

    size_t num_out_map_total;
    size_t num_out_map_preset;
    size_t indel_phasing_neither_flank;
} // namespace global

void set_chr_phase_m()
{
    for (const auto & s : global::chr_phase.get())
    {
        size_t i = s.find(':');
        if (i != s.size() - 2 or (s[s.size() - 1] != '0' and (s[s.size() - 1] != '1')))
        {
            cerr << "error parsing phased chromosome [" << s << "]: format should be <chr:[01]>" << endl;
            abort();
        }
        int phase = s[s.size() - 1] - '0';
        LOG("main", info) << "chromosome [" << s.substr(0, i) << "] set to phase [" << phase << "]" << endl;
        global::chr_phase_m[s.substr(0, i)] = phase;
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
    auto rec_p = bcf_init1();
    int * dat = nullptr; // = new int [bcf_hdr_nsamples(hdr_p)*2];
    int dat_size = 0;
    // restrict processing to given sample_id
    ret = bcf_hdr_set_samples(hdr_p, global::sample_id.get().c_str(), 0);
    if (ret != 0)
    {
        cerr << "bcf_hdr_set_samples() error: status=" << ret << endl;
        exit(1);
    }
    // main loop
    while (bcf_read1(f_p, hdr_p, rec_p) >= 0)
    {
        bcf_unpack(rec_p, BCF_UN_ALL);
        Het_Variation v(hdr_p, rec_p, dat, dat_size, true);
        if (global::chr_phase_m.count(v.chr_name()) > 0 or not v.is_het())
        {
            continue;
        }
        v.load_flanks(global::faidx_p, global::flank_len);
        LOG("variations", debug) << v << endl;
        auto res = global::var_m[v.chr_name()].insert(move(v));
        if (not res.second)
        {
            cerr << "duplicate variation:\n" << v << endl;
            exit(EXIT_FAILURE);
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

/*
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
*/

void close_out_map_files();
void open_out_map_files(const string&);

void implement_decision(const Mapping * m_p, int decision, vector< pair< const Het_Variation *, int > > decision_v)
{
    assert(decision == 0 or decision == 1);
    ++global::num_out_map_total;

    if (not global::decision_fn.get().empty())
    {
        global::decision_ofs
            << m_p->query_name() << '\t'
            << m_p->chr_name() << '\t'
            << decision << '\t';
        bool first = true;
        for (const auto & p : decision_v)
        {
            global::decision_ofs << (not first? ";" : "") << p.first->rf_start() + 1;
            first = false;
        }
        global::decision_ofs << '\t';
        first = true;
        for (const auto & p : decision_v)
        {
            global::decision_ofs << (not first? ";" : "") << p.second;
            first = false;
        }
        global::decision_ofs << endl;
    }
    int res = sam_write1(global::crt_out_map_file_p[decision], global::map_hdr_p, m_p->rec_p());
    if (res < 0)
    {
        cerr << "error in sam_write1(): status=" << res << endl;
        exit(EXIT_FAILURE);
    }
}

int get_unpaired_decision(int decision)
{
    ++global::num_out_map_unpaired_total;
    if (decision < 0)
    {
        if (decision == -1)
        {
            ++global::num_out_map_unpaired_inconclusive;
        }
        else
        {
            ++global::num_out_map_unpaired_random;
        }
        decision = lrand48() % 2;
    }
    else
    {
        ++global::num_out_map_unpaired_phased;
    }
    return decision;
}

int get_paired_decision(int decision, int mp_decision, int frag_decision)
{
    ++global::num_out_frag_total;
    if ((decision < -1) and (mp_decision < -1))
    {
        assert(frag_decision < -1);
        ++global::num_out_frag_hets_neither_side;
    }
    else if ((decision < -1) or (mp_decision < -1))
    {
        // one side has no hets
        assert(frag_decision == (decision >= -1? decision : mp_decision));
        ++global::num_out_frag_hets_single_side;
        if (frag_decision == -1)
        {
            ++global::num_out_frag_hets_single_side_inconclusive;
        }
    }
    else
    {
        // both sides have hets
        ++global::num_out_frag_hets_both_sides;
        if (frag_decision == -1)
        {
            ++global::num_out_frag_hets_both_sides_inconclusive;
        }
    }
    if (frag_decision < 0)
    {
        ++global::num_out_frag_random;
        frag_decision = lrand48() % 2;
    }
    return frag_decision;
}

string get_out_map_fn(const string& chr_name, int phase)
{
    ostringstream os;
    os << global::out_map_fn.get() << "." << chr_name << "." << phase << ".bam";
    return os.str();
}

void close_out_map_files(const string & chr_name)
{
    assert(global::crt_out_map_file_p[0] and global::crt_out_map_file_p[1]);
    // first, flush any paired reads that are currently stored
    for (const auto & p : global::mp_store_m)
    {
        const Mapping * m_p;
        int decision;
        vector< pair< const Het_Variation *, int > > decision_v;
        tie(m_p, decision, decision_v) = p.second;
        if (not m_p)
        {
            ++global::num_out_frag_missing_unmapped;
        }
        else
        {
            ++global::num_out_frag_missing_mapped;
            decision = get_paired_decision(decision, -2, decision);
            implement_decision(m_p, decision, decision_v);
            bam_destroy1(m_p->rec_p());
            delete m_p;
        }
    }
    global::mp_store_m.clear();
    // close files
    hts_close(global::crt_out_map_file_p[0]);
    hts_close(global::crt_out_map_file_p[1]);
    global::crt_out_map_file_p[0] = nullptr;
    global::crt_out_map_file_p[1] = nullptr;
    LOG("main", info) << "closed output files for chr=[" << chr_name << "]" << endl;
}

void open_out_map_files(const string& chr_name)
{
    assert(not global::crt_out_map_file_p[0] and not global::crt_out_map_file_p[1]);
    global::crt_out_map_file_p[0] = hts_open(get_out_map_fn(chr_name, 0).c_str(), "wb");
    global::crt_out_map_file_p[1] = hts_open(get_out_map_fn(chr_name, 1).c_str(), "wb");
    int res = sam_hdr_write(global::crt_out_map_file_p[0], global::map_hdr_p);
    if (res != 0)
    {
        cerr << "error in sam_hdr_write(): status=" << res << endl;
        exit(EXIT_FAILURE);
    }
    res = sam_hdr_write(global::crt_out_map_file_p[1], global::map_hdr_p);
    if (res != 0)
    {
        cerr << "error in sam_hdr_write(): status=" << res << endl;
        exit(EXIT_FAILURE);
    }
    LOG("main", info) << "opened output files for chr=[" << chr_name << "]" << endl;
}

int get_majority_vote(const vector< pair< const Het_Variation *, int > > & decision_v)
{
    int cnt[2] = {0, 0};
    for (const auto & p : decision_v)
    {
        if (p.second >= 0 and p.second < 2)
        {
            // use GQ as factor: put less weight on problematic sites
            int factor = max(1 + p.first->gq(), 1);
            cnt[p.second] += factor;
        }
    }
    return cnt[0] == cnt[1]? -1 : cnt[1] > cnt[0];
}

/**
 * Process one mapping.
 * Returns:
 *   0: if bam record can be reused
 *   1: bam record was saved for deferred decision
 */
int process_mapping(bam1_t * rec_p)
{
    Mapping m(global::map_hdr_p, rec_p);
    bool mp_stored = global::mp_store_m.count(m.query_name()) > 0;
    assert(not mp_stored or m.treat_as_paired());
    //bam1_t * mp_rec_p = nullptr;
    const Mapping * mp_m_p = nullptr;
    int mp_decision = -2;
    vector< pair< const Het_Variation *, int > > mp_decision_v;
    if (mp_stored)
    {
        tie(mp_m_p, mp_decision, mp_decision_v) = global::mp_store_m.at(m.query_name());
    }
    LOG("mappings", debug)
        << m << "\t"
        << (mp_stored? "mp_stored" : "mp_not_stored") << "\t"
        << (mp_stored? (mp_m_p? "non_null_rec" : "null_rec") : "") << "\t"
        << (mp_stored? mp_decision : -5) << "\t"
        << endl;

    // update input mapping counts
    ++global::num_in_map_total;
    if (not m.is_primary())
    {
        // drop non-primary mappings
        ++global::num_in_map_nonprimary;
        return 0;
    }
    if (not m.is_paired())
    {
        assert(not mp_stored);
        ++global::num_in_map_unpaired;
        if (not m.is_mapped())
        {
            // drop unpaired & unmapped
            ++global::num_in_map_unpaired_unmapped;
            return 0;
        }
        else
        {
            ++global::num_in_map_unpaired_mapped;
        }
    }
    else // primary, paired
    {
        if (m.is_mapped())
        {
            ++global::num_in_map_paired_mapped;
        }
        else
        {
            ++global::num_in_map_paired_unmapped;
        }
        if (not m.is_mapped() and not m.mp_is_mapped())
        {
            // drop paired & both unmapped
            ++global::num_in_map_paired_both_unmapped;
            return 0;
        }
        if (m.is_mapped() and m.mp_is_mapped())
        {
            ++global::num_in_map_paired_both_mapped;
            if(m.chr_name() != m.mp_chr_name())
            {
                assert(not mp_stored);
                ++global::num_in_map_paired_both_mapped_diff_chr;
            }
        }
    }
    assert(m.is_mapped() or m.treat_as_paired());

    // compute decision for current mapping
    int decision;
    vector< pair< const Het_Variation *, int > > decision_v;
    if (not m.is_mapped())
    {
        decision = -2; // unmapped
    }
    else
    {
        if (global::chr_phase_m.count(m.chr_name()))
        {
            ++global::num_out_map_preset;
            decision = global::chr_phase_m.at(m.chr_name()); // preset
        }
        else if (not global::var_m.count(m.chr_name()))
        {
            decision = -3; // not spanning any hets
        }
        else
        {
            // compute range of variations spanned by this mapping
            auto it_start = global::var_m.at(m.chr_name()).lower_bound(Het_Variation(m.rf_start()));
            while (it_start != global::var_m.at(m.chr_name()).begin()
                   and prev(it_start)->rf_end() > m.rf_start())
            {
                --it_start;
            }
            auto it_end = global::var_m.at(m.chr_name()).lower_bound(Het_Variation(m.rf_end()));
            if (it_start != it_end)
            {
                // mapping spans at least one het
                for (auto it = it_start; it != it_end; ++it)
                {
                    decision_v.push_back(make_pair(&*it, Phaser::phase(m, *it)));
                }
                decision = get_majority_vote(decision_v);
            }
            else
            {
                decision = -3; // not spanning any hets
            }
        }
    }

    if (not m.treat_as_paired())
    {
        decision = get_unpaired_decision(decision);
        implement_decision(&m, decision, decision_v);
        return 0;
    }
    else
    {
        if (not mp_stored)
        {
            if (m.is_mapped() and not m.mp_is_mapped())
            {
                // paired & 1st & 1st mapped & 2nd unmapped
                // => implement decision & store it
                decision = get_paired_decision(decision, -2, decision);
                implement_decision(&m, decision, decision_v);
                global::mp_store_m[m.query_name()] = make_tuple(nullptr, decision, decision_v);
                return 0;
            }
            else
            {
                // paired & 1st & (1st unmapped or 2nd mapped)
                // => defer decision
                string query_name = m.query_name();
                global::mp_store_m[query_name] = make_tuple(new Mapping(move(m)), decision, decision_v);
                return 1;
            }
        }
        else
        {
            // paired & 2nd
            global::mp_store_m.erase(m.query_name());
            if (mp_m_p)
            {
                // decision for 1st read was deferred
                // merge decision vectors:
                map< const Het_Variation *, int, Het_Variation_Ptr_Comp > decision_m(
                    decision_v.begin(), decision_v.end());
                for (const auto & p : mp_decision_v)
                {
                    if (decision_m.count(p.first) == 0)
                    {
                        decision_m.insert(p);
                    }
                    else if (decision_m[p.first] != p.second)
                    {
                        ++(p.first->frag_conflicting);
                        decision_m.erase(p.first);
                    }
                }
                decltype(decision_v) merged_decision_v(decision_m.begin(), decision_m.end());
                int frag_decision = (not decision_v.empty() or not mp_decision_v.empty()
                                     ? get_majority_vote(merged_decision_v)
                                     : -3);
                frag_decision = get_paired_decision(decision, mp_decision, frag_decision);
                implement_decision(mp_m_p, frag_decision, merged_decision_v);
                implement_decision(&m, frag_decision, merged_decision_v);
                bam_destroy1(mp_m_p->rec_p());
                delete mp_m_p;
                return 0;
            }
            else
            {
                // decision for 1st read was not deferred
                assert(m.is_paired() and m.mp_is_mapped() and not m.is_mapped() and decision_v.empty());
                implement_decision(&m, mp_decision, mp_decision_v);
                return 0;
            }
        }
    }
}

void process_mappings_chromosome(const string & chr_name)
{
    bam1_t * rec_p = bam_init1();
    int tid = bam_name2id(global::map_hdr_p, chr_name.c_str());
    if (tid < 0)
    {
        cerr << "chromosome [" << chr_name << "] not found in BAM header" << endl;
        exit(EXIT_FAILURE);
    }
    auto map_itr_p = sam_itr_queryi(global::map_idx_p, tid, 0, numeric_limits< int >::max());
    auto get_next_record = [&] () { return sam_itr_next(global::map_file_p, map_itr_p, rec_p); };

    open_out_map_files(chr_name);
    while (get_next_record() >= 0)
    {
        int res = process_mapping(rec_p);
        if (res == 1)
        {
            rec_p = bam_init1();
        }
    }
    close_out_map_files(chr_name);

    hts_itr_destroy(map_itr_p);
    bam_destroy1(rec_p);
}

void process_mappings()
{
    global::map_file_p = hts_open(global::map_fn.get().c_str(), "r");
    global::map_hdr_p = sam_hdr_read(global::map_file_p);
    global::map_idx_p = sam_index_load(global::map_file_p, global::map_fn.get().c_str());
    if (not global::map_idx_p)
    {
        cerr << "could not load BAM index for [" << global::map_fn.get() << "]" << endl;
        exit(EXIT_FAILURE);
    }
    if (global::chr.get().empty())
    {
        for (int i = 0; i < global::map_hdr_p->n_targets; ++i)
        {
            process_mappings_chromosome(string(global::map_hdr_p->target_name[i]));
        }
    }
    else
    {
        process_mappings_chromosome(global::chr);
    }
    hts_idx_destroy(global::map_idx_p);
    bam_hdr_destroy(global::map_hdr_p);
    hts_close(global::map_file_p);
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
        << "....hets_neither_side: " << global::num_out_frag_hets_neither_side << endl
        << "....hets_single_side: " << global::num_out_frag_hets_single_side << endl
        << "......hets_single_side_inconclusive: " << global::num_out_frag_hets_single_side_inconclusive << endl
        << "....hets_both_sides: " << global::num_out_frag_hets_both_sides << endl
        << "......hets_both_sides_inconclusive: " << global::num_out_frag_hets_both_sides_inconclusive << endl
        << "....random: " << global::num_out_frag_random << endl
        << "....missing_mapped: " << global::num_out_frag_missing_mapped << endl
        << "....missing_unmapped: " << global::num_out_frag_missing_unmapped << endl

        << "..preset: " << global::num_out_map_preset << endl
        ;
}

void print_var_stats(ostream & os)
{
    if (global::chr.get().empty())
    {
        for (const auto & s : global::var_m)
        {
            for(const auto & v : s.second)
            {
                os << v << endl;
            }
        }
    }
    else
    {
        if (global::var_m.count(global::chr))
        {
            for (const auto & v : global::var_m.at(global::chr))
            {
                os << v << endl;
            }
        }
    }
}

void real_main()
{
    set_chr_phase_m();
    load_reference_index();
    load_variations();
    //print_variations(clog);
    process_mappings();
    print_stats(cout);
    if (not global::var_stats_fn.get().empty())
    {
        strict_fstream::fstream tmp_fs(global::var_stats_fn, ios_base::out);
        print_var_stats(tmp_fs);
    }
    // cleanup
    fai_destroy(global::faidx_p);
}

int main(int argc, char * argv[])
{
    global::cmd_parser.parse(argc, argv);
    // set log levels
    logger::Logger::set_levels_from_options(global::log_level, &clog);
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
    if (not global::decision_fn.get().empty())
    {
        global::decision_ofs.open(global::decision_fn);
    }
    real_main();
    if (not global::decision_fn.get().empty())
    {
        global::decision_ofs.close();
    }
}
