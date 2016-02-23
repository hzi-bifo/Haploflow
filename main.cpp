#include "deBruijnGraph.h"
#include <ctime> //debug
#include <unordered_set> //tmp
#include <algorithm>

int main (int argc, char* argv[])
{
	//clock_t t = clock();
	deBruijnGraph dbg = deBruijnGraph(std::string(argv[1]),atoi(argv[2]));
	//std::cout << (clock() - t)/1000000. << std::endl;
	//t = clock();
	std::vector<std::string> sources = dbg.getSources();
	auto s_tmp = dbg.getSinks();
	std::unordered_set<std::string> sinks(s_tmp.begin(), s_tmp.end());
	std::unordered_map<std::string, std::string> junctions; // maps vertices to their next junction
	for (const auto& source : sources)
	{
		std::string curr = source;
		std::string next = dbg.find_next_junction(curr);
		while (junctions.find(next) == junctions.end() and curr != next)
		{
			junctions.emplace(curr, next);
			curr = next;
			next = dbg.find_next_junction(curr);
		}
	}
	for (const auto& elem : junctions)
	{
		std::string neigh = elem.first.substr(1);
		if (neigh + "A" != elem.second and neigh + "C" != elem.second and neigh + "G" != elem.second and neigh + "T" != elem.second)
			std::cout << elem.first << " --> " << elem.second << std::endl; 
	}
	
	//std::cout << (clock() - t)/1000000. << std::endl;
	//do something
	return 0;
}
