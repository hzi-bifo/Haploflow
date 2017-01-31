#include "deBruijnGraph.h"
#include "UnitigGraph.h"

int main (int argc, char* argv[])
{
	clock_t t = clock();
	std::cerr << "Building deBruijnGraph..." << std::endl;
	deBruijnGraph dbg = deBruijnGraph(std::string(argv[1]), atoi(argv[2]));
	std::cerr << "Building deBruijnGraph took " << (clock() - t)/1000000. << " seconds." << std::endl;
	UnitigGraph ug(dbg);
	return 0;
}
