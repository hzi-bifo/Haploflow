#include "UnitigGraph.h"

UnitigGraph::UnitigGraph(const deBruijnGraph& dbg)
{
	// pair with succ > 1 (first), pred > 1 (second)
	auto& unbalanced = dbg.getUnbalanaced();
	for (const auto& v : unbalanced.first)
	{
		std::stack<Vertex> s;
		for (const auto& n : v.get_successors(false))
		{
			std::string nkmer = v.kmer.substring(1) + n;
			Vertex vnext = dbg.find(nkmer);
		}
		while (s.size() > 0)
		{
			Vertex curr = s.top();
			s.pop();
			std::string edge = "";
			
		}
	}
}
