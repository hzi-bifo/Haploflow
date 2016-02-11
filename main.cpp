#include "deBruijnGraph.h"
#include <ctime> //debug

int main (int argc, char* argv[])
{
	clock_t t = clock();
	deBruijnGraph* dbg = new deBruijnGraph(std::string(argv[1]),atoi(argv[2]));
	//std::cout << (clock() - t)/1000000. << std::endl;
	t = clock();
	const auto& sources = dbg->get_terminals(false); // get sources
	//std::cout << (clock() - t)/1000000. << std::endl;
	t = clock();
	int i = 1;
	for (const auto& v : sources)
	{
		if (!dbg->bfs(v,i++,true))
			i--;
	}
	std::vector<int> count_prev(i - 1);
	for (const auto& v : dbg-> graph_)
		if (v.second[8] > 0)
			count_prev[v.second[8] - 1] += 1;
	for (const auto& u : count_prev)
		std::cout << u << " ";
	std::cout << std::endl;
	for (const auto& v : dbg->graph_)
	{
		if (v.second[8] == 0 and !dbg->bfs(v.first,i++,true))
			i--;
	}
	std::vector<int> count(i - 1);
	for (const auto& v : dbg-> graph_)	
		count[v.second[8] - 1] += 1;
	for (const auto& u : count)
		std::cout << u << " ";
	//std::cout << (clock() - t)/1000000. << std::endl;
	//dbg->printGraph();
	std::cout << i << std::endl;
	delete dbg;
	//do something
	return 0;
}
