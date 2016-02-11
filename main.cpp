#include "deBruijnGraph.h"
#include <ctime> //debug

int main (int argc, char* argv[])
{
	clock_t t = clock();
	deBruijnGraph* dbg = new deBruijnGraph(std::string(argv[1]),atoi(argv[2]));
	std::cout << (clock() - t)/1000000. << std::endl;
	t = clock();
	const auto& list = dbg->get_terminals(false); // get sources
	std::cout << (clock() - t)/1000000. << std::endl;
	delete dbg;
	//do something
	return 0;
}
