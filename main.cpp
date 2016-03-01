#include "deBruijnGraph.h"
#include <ctime> //debug

int main (int argc, char* argv[])
{
	clock_t t = clock();
	deBruijnGraph dbg = deBruijnGraph(std::string(argv[1]),atoi(argv[2]));
	std::cout << (clock() - t)/1000000. << std::endl;
	
	//t = clock();
	//auto c = dbg.find_all_junctions();
	//auto seqs = dbg.getScaffolds(c);
	//for (const auto& seq : seqs)
	//	std::cout << seq << std::endl;

	//std::cout << (clock() - t)/1000000. << std::endl;
	//do something
	return 0;
}
