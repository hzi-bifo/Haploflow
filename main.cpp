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
	for (const auto& v : sources)
	{
		if (!dbg.mark_ccs(v,i++,true))
			i--;
	}
	for (const auto& v : dbg.graph_)
	{
		if (v.second[8] == 0 and !dbg.mark_ccs(v.first,i++,true))
			i--;
	}
	/*debug
	std::cout << dbg.graph_.size() << std::endl;
	std::unordered_map<int,int> ids;
	for (const auto& v : dbg.graph_)
	{
		if (!ids.emplace(v.second[8],1).second)
			ids[v.second[8]]++;
	}
	for (const auto& elem : ids)
		std::cout << elem.first << ": " << elem.second << std::endl;
	end debug*/
		
	std::set<int> id_list;
	std::vector<std::string> sequences;
	int num = 0;
	for (const auto& v : sources)
	{
		if (id_list.find(dbg.graph_[v][8]) == id_list.end())
		{
			dbg.extractSequence(v, num++);
			id_list.insert(dbg.graph_[v][8]);
		}
	}
	for (const auto& v : sinks)
	{
		if (id_list.find(dbg.graph_[v][8]) == id_list.end())
		{
			dbg.extractSequence(v, num++);
			id_list.insert(dbg.graph_[v][8]);
		}
	}
	for (const auto& v : dbg.graph_)
	{
		if (id_list.find(v.second[8]) == id_list.end())
		{
			dbg.extractSequence(v.first, num++);
			id_list.insert(v.second[8]);
		}
	}
	i = 0;
	//std::cout << (clock() - t)/1000000. << std::endl;
	//dbg.printGraph();
	//do something
	return 0;
}
