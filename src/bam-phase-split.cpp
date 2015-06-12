#include <iostream>
#include <cstdlib>
#include <time.h>

#include <htslib/vcf.h>
#include <tclap/CmdLine.h>
#include "logger.hpp"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "missing_version"
#endif

using namespace std;

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
    ValueArg< string > vcf_fn("", "vcf", "VCF file.", false, "", "file", cmd_parser);
    ValueArg< string > bam_fn("", "bam", "BAM file.", false, "", "file", cmd_parser);
} // namespace global

void real_main()
{}

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
