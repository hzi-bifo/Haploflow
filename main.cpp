#include "deBruijnGraph.h"
#include <ctime> //debug
#include <set> //tmp

int main (int argc, char* argv[])
{
	clock_t t = clock();
	deBruijnGraph dbg = deBruijnGraph(std::string(argv[1]),atoi(argv[2]));
	//std::cout << (clock() - t)/1000000. << std::endl;
	t = clock();
	std::vector<std::string> sources = dbg.get_terminals(false); // get sources
	//std::cout << (clock() - t)/1000000. << std::endl;
	t = clock();
	int i = 1;
	for (const auto& v : sources)
	{
		if (!dbg.bfs(v,i++,true))
			i--;
	}
	std::vector<int> count_prev(i - 1);
	for (const auto& v : dbg.graph_)
		if (v.second[8] > 0)
			count_prev[v.second[8] - 1] += 1;
	for (const auto& v : dbg.graph_)
	{
		if (v.second[8] == 0 and !dbg.bfs(v.first,i++,true))
			i--;
	}
		
	std::set<int> id_list;
	std::vector<std::string> sequences;
	for (const auto& v : sources)
	{
		if (id_list.find(dbg.graph_[v][8]) == id_list.end())
		{
			sequences.push_back(dbg.extractSequence(v));
			id_list.insert(dbg.graph_[v][8]);
		}
	}
	for (const auto& s : sequences)
		std::cout << s << std::endl;
	//std::cout << (clock() - t)/1000000. << std::endl;
	//dbg.printGraph();
	//do something
	return 0;
}
