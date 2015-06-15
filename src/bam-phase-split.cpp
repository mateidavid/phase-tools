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
           << v.flank_seq[0] << "\t" << v.flank_seq[1] << endl;
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
    // mappings
    htsFile * map_file_p;
    bam_hdr_t * map_hdr_p;
    string crt_chr;
    htsFile * crt_out_map_file_p[2];
    // mate pair decision
    map< string, int > mp_dec_m;
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

int get_phase(int rf_start, int rf_len, const string& cigar_string, const DNA_Sequence& seq, const Het_Variation& v)
{
    // for now, only deal with SNPs
    if (v.allele_seq_v[0].size() != 1 or v.allele_seq_v[v.gt[0]].size() != 1 or v.allele_seq_v[v.gt[1]].size() != 1)
    {
        return -1;
    }
    //TODO
    return -1;
}

void process_mapping(const bam1_t * rec_p)
{
    string query_name = bam_get_qname(rec_p);
    bool is_paired = rec_p->core.flag & BAM_FPAIRED;
    bool is_mapped = not (rec_p->core.flag & BAM_FUNMAP);
    bool is_mp_mapped = not (rec_p->core.flag & BAM_FMUNMAP);
    bool is_primary = not (rec_p->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY));
    string chr_name = global::map_hdr_p->target_name[rec_p->core.tid];
    int rf_start = rec_p->core.pos;
    auto cigar_ptr = bam_get_cigar(rec_p);
    string cigar_string;
    for (int i = 0; i < rec_p->core.n_cigar; ++i)
    {
        ostringstream oss;
        oss << bam_cigar_oplen(cigar_ptr[i]);
        cigar_string += oss.str() + bam_cigar_opchr(cigar_ptr[i]);
    }
    int rf_len = bam_cigar2rlen(rec_p->core.n_cigar, cigar_ptr);
    auto seq_ptr = bam_get_seq(rec_p);
    DNA_Sequence seq;
    for (int i = 0; i < rec_p->core.l_qseq; ++i)
    {
        seq += seq_nt16_str[bam_seqi(seq_ptr, i)];
    }
    bool mp_seen = global::mp_dec_m.count(query_name) > 0;
    int mp_decision = -1;
    if (mp_seen)
    {
        mp_decision = global::mp_dec_m.at(query_name);
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
        << cigar_string << "\t"
        << seq << endl;

    // compute range of variations spanned by this mapping
    int decision = -2; // unmapped
    if (is_mapped)
    {
        decision = -1; // not spanning any hets
        if (global::var_m.count(chr_name) > 0)
        {
            Het_Variation search_key;
            search_key.rf_start = rf_start;
            auto it_start = global::var_m.at(chr_name).lower_bound(search_key);
            search_key.rf_start = rf_start + rf_len;
            auto it_end = global::var_m.at(chr_name).lower_bound(search_key);
            if (it_start != it_end)
            {
                // mapping spans at least one het
                // for now: compute phasing using only first phased het
                auto it = it_start;
                for (; it != it_end and not it->is_phased; ++it);
                if (it != it_end)
                {
                    int phase = get_phase(rf_start, rf_len, cigar_string, seq, *it);
                    if (phase != -1)
                    {
                        decision = phase;
                    }
                }
            }
        }
    }

    if (not mp_seen)
    {
        global::mp_dec_m[query_name] = decision;
    }
    else
    {
        global::mp_dec_m.erase(query_name);
    }
}

void process_mappings()
{
    global::map_file_p = hts_open(global::map_fn.get().c_str(), "r");
    global::map_hdr_p = sam_hdr_read(global::map_file_p);

    auto rec_p = bam_init1();
    while (sam_read1(global::map_file_p, global::map_hdr_p, rec_p) >= 0)
    {
        process_mapping(rec_p);
    }
    bam_destroy1(rec_p);

    bam_hdr_destroy(global::map_hdr_p);
    hts_close(global::map_file_p);
}

void real_main()
{
    set_chr_phase_m();
    load_reference_index();
    load_variations();
    //print_variations(clog);

    process_mappings();

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
