#include "deBruijnGraph.h"
#include "UnitigGraph.h"
#include "boost/program_options.hpp" 

int main (int argc, char* argv[])
{
    namespace go = boost::program_options;
    go::options_description desc("HaploFlow parameters");
    unsigned int k;
    desc.add_options()
        ("help", "Produce this help message")
        ("read-file, R", go::value<std::string>(), "read file (fastq)")
        ("dump-file, D", go::value<std::string>(), "deBruijn graph dump file produced by HaploFlow")
        ("log", go::value<std::string>(), "log file (default: standard out)")
        ("k, K", go::value<unsigned int>(&k)->default_value(41), "k-mer size, default 41, please use an odd number")
        ("cov, C", go::value<std::string>(), "folder for coverage diagrams, will be created if not present")
        ("graphs, G", go::value<std::string>(), "folder for assembly graphs, will be created if not present")
        ("error-rate, E", go::value<float>(), "percentage filter for erroneous kmers - kmers appearing less than relatively e% will be ignored")
        ("create-dump", go::value<std::string>(), "create dump of the deBruijn graph. WARNING this file may be huge")
    ;
    go::positional_options_description p;
    p.add("read-file", -1);
    go::variables_map vm;
    go::store(go::parse_command_line(argc, argv, desc), vm);
    go::notify(vm);
    if (vm.count("help"))
    {
        std::cout << desc << std::endl;
        return 0;
    }
	clock_t t = clock();
	clock_t t_start = clock();
	std::cerr << "Building deBruijnGraph..." << std::endl;
	//deBruijnGraph *dbg = new deBruijnGraph(std::string(argv[1]), atoi(argv[2])); // read file + k
	deBruijnGraph *dbg = new deBruijnGraph(std::string(argv[1])); // dump file
    std::cerr << "Building deBruijnGraph took " << (clock() - t)/1000000. << " seconds." << std::endl;
    t = clock();
    UnitigGraph ug = UnitigGraph(*dbg, std::string(argv[5]), atof(argv[3])); //argv[2] is k
    delete dbg;
    //t = clock();
    std::cerr << "Assembling..." << std::endl;
    ug.assemble(argv[4], atof(argv[3])); // path to graph files folder
    std::cerr << "Assembly took " << (clock() - t)/1000000. << " seconds" << std::endl;
	std::cerr << "The complete assembly process took " << (clock() - t_start)/1000000. << " seconds." << std::endl;
    return 0;
    /**/
}
