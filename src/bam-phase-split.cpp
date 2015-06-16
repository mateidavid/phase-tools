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

#include "logger.hpp"
#include "RC_Sequence.hpp"
#include "Cigar.hpp"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "missing_version"
#endif

using namespace std;

typedef rc_sequence::Sequence< string > DNA_Sequence;

class Het_Variation
{
public:
    DNA_Sequence flank_seq[2];
    vector< DNA_Sequence > allele_seq_v;
    int rf_start;
    int rf_len;
    int gt[2];
    bool is_phased;
    bool is_snp; // true iff all of (potentially 3) alleles: ref, gt[0], gt[1] are length 1
    unsigned long fragment_count[2];

    friend bool operator < (const Het_Variation& lhs, const Het_Variation& rhs)
    {
        return lhs.rf_start < rhs.rf_start;
    }

    friend ostream& operator << (ostream& os, const Het_Variation& v)
    {
        os << v.rf_start + 1 << "\t";
        for (size_t i = 0; i < v.allele_seq_v.size(); ++i)
        {
            os << v.allele_seq_v[i] << ((i > 0 and i < v.allele_seq_v.size()-1)? "," : "\t");
        }
        os << v.gt[0] << (v.is_phased? "|" : "/") << v.gt[1] << "\t"
           << v.flank_seq[0] << "\t" << v.flank_seq[1] << "\t"
           << (v.is_snp? "snp" : "indel") << endl;
        return os;
    }
}; // class Variation

namespace global
{
    using namespace TCLAP;

    string description =
        "Given a BAM file and a partly phased VCF file, split the mappings in the BAM file "
        "into several other BAM files. Specifically, for every phase set (PS) defined in the VCF file, "
        "and for every autosome (1-22), there will be 2 BAM output files (1 per haplotype). "
        "In addition, there will be 2 more overflow BAM output files for reads not mapped to autosomes.";
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
    ValueArg< string > map_fn("", "map", "Mappings (BAM) file.", true, "", "file", cmd_parser);
    ValueArg< string > out_map_fn("", "out", "Output file prefix. Files names will be <prefix>.<chr>.[01].bam", true, "", "file", cmd_parser);
    ValueArg< string > sample_id("", "sample", "Sample Id.", true, "", "string", cmd_parser);
    //
    // other parameters
    //
    MultiArg< string > chr_phase("", "chr-phase", "Assign all mappings to this chr to given phase.", false, "chr:[01]", cmd_parser);
    ValueArg< int > flank_len("", "flank_len", "Flank length [10].", false, 10, "int", cmd_parser);

    // reference
    faidx_t * faidx_p;
    // phased chromosomes
    map< string, int > chr_phase_m;
    // variations
    map< string, set< Het_Variation > > var_m;
    // input mappings
    htsFile * map_file_p;
    bam_hdr_t * map_hdr_p;
    // output mappings
    string crt_chr;
    htsFile * crt_out_map_file_p[2];
    // mate pair decision
    map< string, pair< bam1_t *, int > > mp_store_m;

    // counts
    size_t num_frag_random_decision;
    size_t num_frag_inconclusive_decision;
    size_t num_frag_single_decision;
    size_t num_frag_concordant_decision;
    size_t num_frag_conflicting_decision;
    size_t num_frag_unmapped;
    size_t num_frag_diff_chr;
    size_t num_map_nonprimary;
    size_t num_map_output;
    size_t num_map_inconclusive_phasing;
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
    int * dat = nullptr;// = new int [bcf_hdr_nsamples(hdr_p)*2];
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
        Het_Variation v;
        bcf_unpack(rec_p, BCF_UN_ALL);
        int n_gt = bcf_get_genotypes(hdr_p, rec_p, &dat, &dat_size);
        string chr_name = bcf_hdr_id2name(hdr_p, rec_p->rid);
        v.rf_start = rec_p->pos;
        v.rf_len = rec_p->rlen;

        /*
        clog << chr_name << "\t" << v.rf_start << "\t" << v.rf_len << "\t";
        clog << "n_allele=" << rec_p->n_allele << "\t";
        for (int i = 0; i < rec_p->n_allele; ++i)
            clog << rec_p->d.allele[i] << ",";
        clog << "\t";

        clog << "dat_size=" << dat_size << "\t" << dat[0] << ":" << dat[1] << "\t";
        if (bcf_gt_is_missing(dat[0]))
            clog << ".";
        else
            clog << bcf_gt_allele(dat[0]);
        if (bcf_gt_is_phased(dat[1]))
            clog << "|";
        else
            clog << "/";
        if (bcf_gt_is_missing(dat[1]))
            clog << ".";
        else
            clog << bcf_gt_allele(dat[1]);
        clog << endl;
        */

        if (global::chr_phase_m.count(chr_name) > 0
            or n_gt != 2
            or bcf_gt_is_missing(dat[0])
            or bcf_gt_is_missing(dat[1])
            //or not bcf_gt_is_phased(dat[1])
            or bcf_gt_allele(dat[0]) == bcf_gt_allele(dat[1]))
        {
            continue;
        }

        {
            char * seq;
            int seq_size;
            seq = faidx_fetch_seq(global::faidx_p, chr_name.c_str(),
                                  v.rf_start - global::flank_len, v.rf_start - 1,
                                  &seq_size);
            v.flank_seq[0] = seq;
            free(seq);
            seq = faidx_fetch_seq(global::faidx_p, chr_name.c_str(),
                                  v.rf_start + v.rf_len, v.rf_start + v.rf_len + global::flank_len - 1,
                                  &seq_size);
            v.flank_seq[1] = seq;
            free(seq);
        }
        for (int i = 0; i < rec_p->n_allele; ++i)
        {
            v.allele_seq_v.emplace_back(rec_p->d.allele[i]);
        }
        for (int i = 0; i < 2; ++i)
        {
            v.gt[i] = bcf_gt_allele(dat[i]);
        }
        v.is_phased = bcf_gt_is_phased(dat[1]);
        v.is_snp = (v.allele_seq_v[0].size() == 1
                    and v.allele_seq_v[v.gt[0]].size() == 1
                    and v.allele_seq_v[v.gt[1]].size() == 1);
        if (not v.is_phased)
        {
            // assign a random phase
            int flip_phase = lrand48() % 2;
            if (flip_phase) swap(v.gt[0], v.gt[1]);
        }
        LOG("variations", debug) << chr_name << "\t" << v;
        auto res = global::var_m[chr_name].insert(move(v));
        if (not res.second)
        {
            cerr << "duplicate variation:\n" << chr_name << "\t" << v;
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

int get_phase_snp(int rf_start, int rf_len, const Cigar& cigar, const DNA_Sequence& seq, const Het_Variation& v)
{
    assert(rf_start <= v.rf_start and v.rf_start + v.rf_len <= rf_start + rf_len);
    // compute the mapped range, on query, of the single reference position
    auto rg = cigar.mapped_range(make_pair(v.rf_start - rf_start, v.rf_start - rf_start + 1), true);
    if (rg.second != rg.first + 1)
    {
        // the changed reference base is mapped to something other than a match; we give up
        return -1;
    }
    else
    {
        // the changed reference base is mapped to seq[rg.first]
        assert(0 <= rg.first and rg.second <= static_cast< int >(seq.size()));
        if (seq[rg.first] == v.allele_seq_v[v.gt[0]][0])
            return 0;
        else if (seq[rg.first] == v.allele_seq_v[v.gt[1]][0])
            return 1;
        else
            return -1;
    }
}

int get_phase_indel(int, int, const Cigar&, const DNA_Sequence&, const Het_Variation&)
{
    return -1;
}

int get_phase(int rf_start, int rf_len, const Cigar& cigar, const DNA_Sequence& seq, const Het_Variation& v)
{
    // for now, only deal with SNPs
    return v.is_snp
        ? get_phase_snp(rf_start, rf_len, cigar, seq, v)
        : get_phase_indel(rf_start, rf_len, cigar, seq, v);
}

string get_out_map_fn(const string& chr_name, int phase)
{
    ostringstream os;
    os << global::out_map_fn.get() << "." << chr_name << "." << phase << ".bam";
    return os.str();
}

void close_out_map_files()
{
    assert((global::crt_out_map_file_p[0] != nullptr) == (global::crt_out_map_file_p[1] != nullptr));
    if (global::crt_out_map_file_p[0])
    {
        hts_close(global::crt_out_map_file_p[0]);
        hts_close(global::crt_out_map_file_p[1]);
    }
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
    global::crt_chr = chr_name;
    LOG("main", info) << "opened output files for chr=[" << chr_name << "]" << endl;
}

void implement_decision(const bam1_t * rec_p, const string& chr_name, int decision)
{
    assert(decision == 0 or decision == 1);
    ++global::num_map_output;

    if (global::crt_chr != chr_name)
    {
        // new chromosome
        close_out_map_files();
        open_out_map_files(chr_name);
    }
    int res = sam_write1(global::crt_out_map_file_p[decision], global::map_hdr_p, rec_p);
    if (res < 0)
    {
        cerr << "error in sam_write1(): status=" << res << endl;
        exit(EXIT_FAILURE);
    }
}

void set_single_decision(int & decision)
{
    if (decision < 0)
    {
        if (decision == -1)
        {
            ++global::num_frag_inconclusive_decision;
        }
        else
        {
            ++global::num_frag_random_decision;
        }
        decision = lrand48() % 2;
    }
}

/**
 * Process one mapping.
 * Returns:
 *   0: if bam record can be reused
 *   1: bam record was saved for deferred decision
 */
int process_mapping(bam1_t * rec_p)
{
    string query_name = bam_get_qname(rec_p);
    bool is_paired = rec_p->core.flag & BAM_FPAIRED;
    bool is_mapped = not (rec_p->core.flag & BAM_FUNMAP);
    bool is_mp_mapped = not (rec_p->core.flag & BAM_FMUNMAP);
    bool is_primary = not (rec_p->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY));
    string chr_name = global::map_hdr_p->target_name[rec_p->core.tid];
    string mp_chr_name;
    if (is_mp_mapped)
    {
        mp_chr_name = global::map_hdr_p->target_name[rec_p->core.mtid];
    }
    int rf_start = rec_p->core.pos;
    Cigar cigar(bam_get_cigar(rec_p), rec_p->core.n_cigar);
    int rf_len = cigar.rf_len();
    auto seq_ptr = bam_get_seq(rec_p);
    DNA_Sequence seq;
    for (int i = 0; i < rec_p->core.l_qseq; ++i)
    {
        seq += seq_nt16_str[bam_seqi(seq_ptr, i)];
    }
    bool mp_seen = global::mp_store_m.count(query_name) > 0;
    int mp_decision = -1;
    if (mp_seen)
    {
        mp_decision = global::mp_store_m.at(query_name).second;
    }

    LOG("mappings", debug)
        << query_name << "\t"
        << (is_paired? "paired" : "unpaired") << "\t"
        << (is_mapped? "mapped" : "unmapped") << "\t"
        << (is_mp_mapped? "mp_mapped" : "mp_unmapped") << "\t"
        << (is_primary? "primary" : "non-primary") << "\t"
        << (mp_seen? "second" : "first") << "\t"
        << mp_decision << "\t"
        << chr_name << "\t"
        << rf_start << "\t"
        << rf_len << "\t"
        << cigar.to_string() << "\t"
        << seq << endl;

    // drop non-primary mappings
    if (not is_primary)
    {
        ++global::num_map_nonprimary;
        return 0;
    }

    // compute range of variations spanned by this mapping
    int decision = -2; // unmapped
    if (is_mapped)
    {
        decision = -3; // not spanning any hets
        if (global::chr_phase_m.count(chr_name))
        {
            decision = global::chr_phase_m.at(chr_name);
        }
        else if (global::var_m.count(chr_name) > 0)
        {
            Het_Variation search_key;
            search_key.rf_start = rf_start;
            auto it_start = global::var_m.at(chr_name).lower_bound(search_key);
            search_key.rf_start = rf_start + rf_len;
            auto it_end = global::var_m.at(chr_name).lower_bound(search_key);
            if (it_start != it_end)
            {
                // mapping spans at least one het
                // for now: compute phasing using only first het
                decision = get_phase(rf_start, rf_len, cigar, seq, *it_start);
                if (decision == -1)
                {
                    ++global::num_map_inconclusive_phasing;
                }
            }
        }
    }

    if (not mp_seen)
    {
        if (not is_mapped)
        {
            if (is_mp_mapped)
            {
                // 1st read
                // 1st read unmapped
                // 2nd read mapped
                // => defer decision
                global::mp_store_m[query_name] = make_pair(rec_p, -2);
                return 1;
            }
            else
            {
                // 1st read
                // 1st read unmapped
                // 2nd read unmapped
                // => dump fragment
                ++global::num_frag_unmapped;
                return 0;
            }
        }
        else
        {
            if (is_mp_mapped and mp_chr_name == chr_name)
            {
                // 1st read
                // 1st mapped
                // 2nd mapped to same chr
                // => defer decision
                global::mp_store_m[query_name] = make_pair(rec_p, decision);
                return 1;
            }
            else
            {
                // 1st read
                // 1st mapped
                // 2nd unmapped or 2nd mapped to different chr
                // => take decision and save it for logging
                if (is_mp_mapped)
                {
                    ++global::num_frag_diff_chr;
                }
                set_single_decision(decision);
                implement_decision(rec_p, chr_name, decision);
                global::mp_store_m[query_name] = make_pair(nullptr, decision);
                return 0;
            }
        }
    }
    else
    {
        if (not is_mapped and not is_mp_mapped)
        {
            // 2nd read
            // 1st & 2nd unmapped
            // => dump fragment
            return 0;
        }
        else
        {
            // 2nd read
            // 1st or 2nd mapped
            // there must be a record of the decision for 1st
            bam1_t * mp_rec_p;
            int mp_decision;
            tie(mp_rec_p, mp_decision) = global::mp_store_m.at(query_name);
            global::mp_store_m.erase(query_name);
            if (mp_rec_p)
            {
                // decision for 1st read was deferred
                assert(not is_mp_mapped or chr_name == mp_chr_name);
                int frag_decision;
                if (decision >= 0 and mp_decision >= 0 and decision != mp_decision)
                {
                    // conflicting decisions
                    ++global::num_frag_conflicting_decision;
                    frag_decision = lrand48() % 2;
                }
                else if (decision < 0 or mp_decision < 0)
                {
                    if (decision < 0 and mp_decision < 0)
                    {
                        if (decision == -1 or mp_decision == -1)
                        {
                            ++global::num_frag_inconclusive_decision;
                        }
                        else
                        {
                            ++global::num_frag_random_decision;
                        }
                        frag_decision = lrand48() % 2;
                    }
                    else
                    {
                        ++global::num_frag_single_decision;
                        frag_decision = decision >= 0? decision : mp_decision;
                    }
                }
                else
                {
                    assert(decision >= 0 and mp_decision >= 0 and decision == mp_decision);
                    ++global::num_frag_concordant_decision;
                    frag_decision = decision;
                }
                implement_decision(mp_rec_p, chr_name, frag_decision);
                implement_decision(rec_p, chr_name, frag_decision);
                bam_destroy1(mp_rec_p);
                return 0;
            }
            else
            {
                // decision for 1st read was not deferred
                assert(is_mp_mapped);
                if (is_mapped)
                {
                    assert(chr_name != mp_chr_name);
                    set_single_decision(decision);
                    implement_decision(rec_p, chr_name, decision);
                }
                else
                {
                    assert(chr_name == mp_chr_name);
                    implement_decision(rec_p, chr_name, mp_decision);
                }
                return 0;
            }
        }
    }
}

void process_mappings()
{
    global::map_file_p = hts_open(global::map_fn.get().c_str(), "r");
    global::map_hdr_p = sam_hdr_read(global::map_file_p);

    bam1_t * rec_p = bam_init1();
    while (sam_read1(global::map_file_p, global::map_hdr_p, rec_p) >= 0)
    {
        int res = process_mapping(rec_p);
        if (res == 1)
        {
            rec_p = bam_init1();
        }
    }
    bam_destroy1(rec_p);

    close_out_map_files();
    bam_hdr_destroy(global::map_hdr_p);
    hts_close(global::map_file_p);
}

void print_stats(ostream & os)
{
    os
        << "fragments_random_decision: " << global::num_frag_random_decision << endl
        << "fragments_inconclusive_decision: " << global::num_frag_inconclusive_decision << endl
        << "fragments_single_decision: " << global::num_frag_single_decision << endl
        << "fragments_concordant_decision: " << global::num_frag_concordant_decision << endl
        << "fragments_conflicting_decision: " << global::num_frag_conflicting_decision << endl
        << "fragments_unmapped: " << global::num_frag_unmapped << endl
        << "fragments_diff_chr: " << global::num_frag_diff_chr << endl
        << "mappings_nonprimary: " << global::num_map_nonprimary << endl
        << "mappings_output: " << global::num_map_output << endl
        << "mappings_inconclusive: " << global::num_map_inconclusive_phasing << endl
        ;
}

void real_main()
{
    set_chr_phase_m();
    load_reference_index();
    load_variations();
    //print_variations(clog);
    process_mappings();
    print_stats(cout);
    // cleanup
    fai_destroy(global::faidx_p);
}

int main(int argc, char * argv[])
{
    global::cmd_parser.parse(argc, argv);
    // set log levels
    Logger::set_levels_from_options(global::log_level, &clog);
    // set random seed
    bool random_seed = false;
    if (global::seed < 0)
    {
        global::seed.get() = time(nullptr);
        random_seed = true;
    }
    srand48(global::seed);
    // print options
    LOG("main", info) << "program: " << global::cmd_parser.getProgramName() << endl;
    LOG("main", info) << "version: " << global::cmd_parser.getVersion() << endl;
    LOG("main", info) << "args: " << global::cmd_parser.getOrigArgv() << endl;
    if (random_seed) LOG("main", info) << "seed: " << global::seed << endl;
    // real main
    real_main();
}
