#include "deBruijnGraph.h"

int main (int argc, char* argv[])
{
	clock_t t = clock();
	deBruijnGraph dbg = deBruijnGraph(std::string(argv[1]), true, atoi(argv[3]));
	//dbg.add_sequence(std::string(argv[2]));
	std::cerr << (clock() - t)/1000000. << std::endl;
	dbg.debug();

	return 0;
}
