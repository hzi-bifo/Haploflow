#include "deBruijnGraph.h"

int main (int argc, char* argv[])
{
	deBruijnGraph dbg = deBruijnGraph(std::string(argv[1]),atoi(argv[2]));
	dbg.debug();

	return 0;
}
