#include "deBruijnGraph.h"
#include "UnitigGraph.h"

int main (int argc, char* argv[])
{
	clock_t t = clock();
	clock_t t_start = clock();
	std::cerr << "Building deBruijnGraph..." << std::endl;
	deBruijnGraph *dbg = new deBruijnGraph(std::string(argv[1]), atoi(argv[2])); // read file + k
	//deBruijnGraph *dbg = new deBruijnGraph(std::string(argv[1])); // dump file
    std::cerr << "Building deBruijnGraph took " << (clock() - t)/1000000. << " seconds." << std::endl;
    t = clock();
    UnitigGraph ug = UnitigGraph(*dbg, std::string(argv[5]), atof(argv[3])); //argv[2] is k
    delete dbg;
    //t = clock();
    std::cerr << "Assembling..." << std::endl;
    ug.assemble(argv[4]); // path to graph files folder
    std::cerr << "Assembly took " << (clock() - t)/1000000. << " seconds" << std::endl;
	std::cerr << "The complete assembly process took " << (clock() - t_start)/1000000. << " seconds." << std::endl;
    return 0;
}
