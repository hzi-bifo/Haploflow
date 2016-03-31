#include "deBruijnGraph.h"
#include <ctime> //debug

int main (int argc, char* argv[])
{
	clock_t t = clock();
	deBruijnGraph dbg = deBruijnGraph(std::string(argv[1]),atoi(argv[2]));
	std::cerr << dbg.getSize() << std::endl;
	std::cerr << (clock() - t)/1000000. << std::endl;
	
	t = clock();
	auto c = dbg.find_all_junctions();
	std::cerr << (clock() - t)/1000000. << std::endl;
	t = clock();
	for (const auto& p : c)
	{
		dbg.getSequences(p.first, p.second);
	}
	std::cout << (clock() - t)/1000000. << std::endl;
	//do something
	return 0;
}
