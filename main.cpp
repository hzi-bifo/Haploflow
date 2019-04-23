#include "deBruijnGraph.h"
#include "UnitigGraph.h"

int main (int argc, char* argv[])
{
	clock_t t = clock();
	std::cerr << "Building deBruijnGraph..." << std::endl;
	//deBruijnGraph dbg = deBruijnGraph(std::string(argv[1]), atoi(argv[2])); // read file + k
	deBruijnGraph *dbg = new deBruijnGraph(std::string(argv[1])); // dump file
    std::cerr << "Building deBruijnGraph took " << (clock() - t)/1000000. << " seconds." << std::endl;
    t = clock();
    std::cerr << "Getting connected components" << std::endl;
    auto dbgs = dbg->split_ccs();
    delete dbg; // we split it up and can delete the original
    std::cerr << "Getting CCs took " << (clock() - t)/1000000. << " seconds" << std::endl;
    //dbg.debug();
    //return 0;
	//std::cout << dbg;
    //return 0;
    //t = clock();
	//dbg.markCycles();
	//std::cerr << "Cycle detection took " << (clock() - t)/1000000. << " seconds." << std::endl;
	//dbg.debug();
    for (auto&& d : dbgs)
    {
	    UnitigGraph ug(*d, atof(argv[2]));
        delete d;
        //t = clock();
        //ug.debug();
        ug.assemble(argv[3]);
    }
	std::cerr << "The complete assembly process took " << (clock() - t)/1000000. << " seconds." << std::endl;
    return 0;
}
