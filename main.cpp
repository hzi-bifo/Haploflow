#include "deBruijnGraph.h"
#include "UnitigGraph.h"

int main (int argc, char* argv[])
{
	clock_t t = clock();
	std::cerr << "Building deBruijnGraph..." << std::endl;
	deBruijnGraph dbg = deBruijnGraph(std::string(argv[1]), atoi(argv[2])); // read file + k
	//deBruijnGraph dbg = deBruijnGraph(std::string(argv[1])); // dump file
    std::cerr << "Building deBruijnGraph took " << (clock() - t)/1000000. << " seconds." << std::endl;
	//std::cout << dbg;
    //return 0;
    t = clock();
	//dbg.markCycles();
	//std::cerr << "Cycle detection took " << (clock() - t)/1000000. << " seconds." << std::endl;
	//dbg.debug();
	UnitigGraph ug(dbg, atof(argv[2]));
    t = clock();
	//ug.debug();
    ug.assemble();
	//std::cerr << "Calculating flow took " << (clock() - t)/1000000. << " seconds." << std::endl;
	return 0;
}
