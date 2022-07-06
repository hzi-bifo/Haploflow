#include "deBruijnGraph.h"
#include "UnitigGraph.h"
#include "boost/program_options.hpp" 
#include "boost/filesystem.hpp"
#include <sys/stat.h>

int main (int argc, char* argv[])
{
    namespace go = boost::program_options;
    go::options_description desc("HaploFlow parameters");
    std::string reads;
    std::string d;
    std::string log;
    unsigned int k;
    std::vector<unsigned int> k_list;
    std::string o;
    float e;
    std::string create_dump;
    std::string from_dump;
    unsigned int strict;
    bool two_strain;
    unsigned int filter;
    int thresh;
    bool l;
    bool true_flow;
    bool debug;
    desc.add_options()
        ("help", "Produce this help message")
        ("read-file, r", go::value<std::string>(&reads)->required(), "read file (fastq)")
        ("dump-file, d", go::value<std::string>(&d), "deBruijn graph dump file produced by HaploFlow")
        ("log", go::value<std::string>(&log), "log file (default: standard out)")
        ("k", go::value<unsigned int>(&k), "Single k-mer mode value, please use an odd number")
        ("k-list, K", go::value<std::vector<unsigned int> >(&k_list)->multitoken()->default_value(std::vector<unsigned int>{41,61,81,101},"41 61 81 101"), "List of values for k, default 41 61 81 101")
        ("out, o", go::value<std::string>(&o)->required(), "folder for output, will be created if not present. WARNING: Old results will get overwritten")
        ("error-rate, e", go::value<float>(&e)->default_value(0.02), "percentage filter for erroneous kmers - kmers appearing less than relatively e% will be ignored")
        ("create-dump", go::value<std::string>(&create_dump), "create dump of the deBruijn graph. WARNING: This file may be huge")
        ("from-dump", go::value<std::string>(&from_dump), "run from a Haploflow dump of the deBruijn graph.")
        ("two-strain, 2", go::value<bool>(&two_strain)->default_value(false), "mode for known two-strain mixtures")
        ("strict, s", go::value<unsigned int>(&strict)->default_value(5), "more strict error correction, should be set to 5 in first run on new data set to reduce run time. Set to 0 if low abundant strains are expected to be present")
        ("filter, f", go::value<unsigned int>(&filter)->default_value(500), "filter contigs shorter than value")
        ("thresh, t", go::value<int>(&thresh)->default_value(-1), "Provide a custom threshold for complex/bad data")
        ("long, l", go::value<bool>(&l)->default_value(false), "Try to maximise contig lengths (might introduce errors)")
        ("true-flow, tf", go::value<bool>(&true_flow)->default_value(false), "Do not perform flow correction, assume perfect flows")
        ("debug, d", go::bool_switch(&debug)->default_value(false), "Report all temporary graphs and coverage histograms")
    ;
    go::positional_options_description p;
    p.add("read-file", -1);
    go::variables_map vm;
    go::store(go::parse_command_line(argc, argv, desc), vm);
    go::notify(vm);
    if (vm.count("help") || argc == 1)
    {
        std::cout << desc << std::endl;
        return 0;
    }
    std::string out = o;
    // UNIX only
    std::string cov = o + "/Coverages/";
    std::string g = o + "/Graphs/";
    std::string contigs = o + "/contigs.fa";
    if (!boost::filesystem::exists(out))
        boost::filesystem::create_directory(o);
    if (!boost::filesystem::exists(g) and debug)
        boost::filesystem::create_directory(g);
    if (!boost::filesystem::exists(cov) and debug)
        boost::filesystem::create_directory(cov);
    if (!debug)
        cov = o;

	clock_t t = clock();
	clock_t t_start = clock();
	std::ofstream logfile;
    if (vm.count("log"))
        logfile.open(log);
    if (debug)
    {
        if (vm.count("log"))
        {
            logfile << "Options used: " << std::endl;
            logfile << "strict " << strict << ", error-rate " << e;
            logfile << ", two-strain " << (two_strain ? "True" : "False") << ", long contigs: " << (l ? "True" : "False");
            logfile << ", filter " << filter << ", threshold " << thresh << std::endl;
            if (vm.count("k"))
            {
                logfile << "Using single value for k: " << k << std::endl;
            }
            else
            {
                logfile << "Using the following values for k: ";
                for (const auto& elem : k_list)
                {
                    logfile << elem << " ";
                }
                logfile << std::endl;
            }
        }
        else
        {
            std::cout << "Options used: " << std::endl;
            std::cout << "strict " << strict << ", error-rate " << e;
            std::cout << ", two-strain " << (two_strain ? "True" : "False") << ", long contigs: " << (l ? "True" : "False");
            std::cout << ", filter " << filter << ", threshold " << thresh << std::endl;
            if (vm.count("k"))
            {
                std::cout << "Using single value for k: " << k << std::endl;
            }
            else
            {
                std::cout << "Using the following values for k: ";
                for (const auto& elem : k_list)
                {
                    std::cout << elem << " ";
                }
                std::cout << std::endl;
            }
        }
    }
    logfile << "Building deBruijnGraph..." << std::endl;
    logfile.close();
    deBruijnGraph* dbg;
	if (vm.count("from-dump"))
	{
        dbg = new deBruijnGraph(from_dump); // dump file
    }
    else
    {
        dbg = new deBruijnGraph(reads, k); // read file + k
    }
    logfile.open(log, std::ofstream::out | std::ofstream::app);
    logfile << "Building deBruijnGraph took " << (clock() - t)/1000000. << " seconds." << std::endl;
    logfile.close();
    t = clock();
    UnitigGraph ug = UnitigGraph(*dbg, cov, log, e, strict, filter, thresh, l, debug); //argv[2] is k
    delete dbg;
    if (true_flow)
    {
        ug.set_debug(); // remove flow correction
    }
    //t = clock();
    logfile.open(log, std::ofstream::out | std::ofstream::app);
    logfile << "Assembling..." << std::endl;
    logfile.close();
    ug.assemble(g, e, contigs, two_strain); // path to graph files folder
    logfile.open(log, std::ofstream::out | std::ofstream::app);
    logfile << "Assembly took " << (clock() - t)/1000000. << " seconds" << std::endl;
	logfile << "The complete assembly process took " << (clock() - t_start)/1000000. << " seconds." << std::endl;
    logfile.close();
    return 0;
    /**/
}
