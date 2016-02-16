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
	std::vector<std::string> sinks = dbg.get_terminals(true); // get sources
	//std::cout << (clock() - t)/1000000. << std::endl;
	t = clock();
	int i = 1;
	for (const auto& v : dbg.graph_)
		if (v.second.cc == 0)
			dbg.mark_ccs(v.first,i++);

	/*debug*/
	std::cout << dbg.graph_.size() << std::endl;
	std::unordered_map<int,int> ids;
	for (const auto& v : dbg.graph_)
	{
		if (!ids.emplace(v.second.cc,1).second)
			ids[v.second.cc]++;
	}
	//for (const auto& elem : ids)
	//	std::cout << elem.first << ": " << elem.second << std::endl;
	/*end debug*/

	//std::cout << (clock() - t)/1000000. << std::endl;
	dbg.printGraph();
	//do something
	return 0;
}
