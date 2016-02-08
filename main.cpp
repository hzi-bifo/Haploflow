#include "deBruijnGraph.h"

int main (int argc, char* argv[])
{
	deBruijnGraph* dbg = new deBruijnGraph(std::string(argv[1]),atoi(argv[2]));
	dbg->printGraph();
	delete dbg;
	return 0;
}
