#include "deBruijnGraph.h"

deBruijnGraph::deBruijnGraph(int k)
{
	k_ = k;
}

deBruijnGraph::deBruijnGraph(std::string filename, int k)
{
	k_ = k;
	std::ifstream infile(filename);
}
