#include <iostream>
#include <cstdlib>
#include <time.h>

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
    map< string, list< Het_Variation > > var_m;
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
        LOG("main", debug) << chr_name << "\t" << v;
        global::var_m[chr_name].push_back(move(v));
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

void process_mappings()
{
    auto map_file_p = hts_open(global::map_fn.get().c_str(), "r");
    auto map_hdr_p = sam_hdr_read(map_file_p);
    auto rec_p = bam_init1();
    while (sam_read1(map_file_p, map_hdr_p, rec_p) >= 0)
    {
        string query_name = bam_get_qname(rec_p);
        bool is_paired = rec_p->core.flag & BAM_FPAIRED;
        bool is_primary = not (rec_p->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY));
        string chr_name = map_hdr_p->target_name[rec_p->core.tid];
        int rf_start = rec_p->core.pos;
        auto cigar_ptr = bam_get_cigar(rec_p);
        string cigar_string;
        for (int i = 0; i < rec_p->core.n_cigar; ++i)
        {
            ostringstream oss;
            oss << bam_cigar_oplen(cigar_ptr[i]);
            cigar_string += oss.str() + bam_cigar_opchr(cigar_ptr[i]);
        }
        auto seq_ptr = bam_get_seq(rec_p);
        DNA_Sequence seq;
        for (int i = 0; i < rec_p->core.l_qseq; ++i)
        {
            switch (bam_seqi(seq_ptr, i))
            {
            case 1:
                seq += 'A';
                break;
            case 2:
                seq += 'C';
                break;
            case 4:
                seq += 'G';
                break;
            case 8:
                seq += 'T';
                break;
            default:
                seq += 'N';
                break;
            }
        }

        clog << query_name << "\t"
             << (is_paired? "paired" : "unpaired") << "\t"
             << (is_primary? "primary" : "non-primary") << "\t"
             << chr_name << "\t"
             << rf_start << "\t"
             << cigar_string << "\t"
             << seq << endl;

    }
    bam_destroy1(rec_p);
    bam_hdr_destroy(map_hdr_p);
    hts_close(map_file_p);
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
