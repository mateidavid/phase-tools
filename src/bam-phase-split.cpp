#include <iostream>
#include <cstdlib>
#include <time.h>

#include <htslib/vcf.h>
#include <tclap/CmdLine.h>
#include "version.h"

#include "logger.hpp"
#include "RC_Sequence.hpp"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "missing_version"
#endif

using namespace std;

class Het_Variation
{
public:
    vector< rc_sequence::Sequence< string > > seq_v;
    int rf_start;
    int rf_len;
    int gt[2];
    bool is_phased;
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
    ValueArg< string > vcf_fn("", "vcf", "VCF file.", true, "", "file", cmd_parser);
    ValueArg< string > bam_fn("", "bam", "BAM file.", true, "", "file", cmd_parser);
    ValueArg< string > sample_id("", "sample", "Sample Id.", true, "", "string", cmd_parser);

    map< string, list< Het_Variation > > var_map;
} // namespace global

void load_variations()
{
    int ret;
    // open vcf file, initialize structs
    auto f_p = bcf_open(global::vcf_fn.get().c_str(), "r");
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

        if (n_gt != 2
            or bcf_gt_is_missing(dat[0])
            or bcf_gt_is_missing(dat[1])
            //or not bcf_gt_is_phased(dat[1])
            or bcf_gt_allele(dat[0]) == bcf_gt_allele(dat[1]))
        {
            continue;
        }

        for (int i = 0; i < rec_p->n_allele; ++i)
        {
            v.seq_v.emplace_back(rec_p->d.allele[i]);
        }
        for (int i = 0; i < 2; ++i)
        {
            v.gt[i] = bcf_gt_allele(dat[i]);
        }
        v.is_phased = bcf_gt_is_phased(dat[1]);
        global::var_map[chr_name].push_back(move(v));
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
}

void real_main()
{
    load_variations();
    for (const auto & p : global::var_map)
    {
        for (const auto & v : p.second)
        {
            clog << p.first << "\t" << v.rf_start + 1 << "\t";
            for (size_t i = 0; i < v.seq_v.size(); ++i)
            {
                clog << v.seq_v[i] << ((i > 0 and i < v.seq_v.size()-1)? "," : "\t");
            }
            clog << v.gt[0] << (v.is_phased? "|" : "/") << v.gt[1] << endl;
        }
    }
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
