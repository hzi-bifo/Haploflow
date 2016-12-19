#include "UnitigGraph.h"

UnitigGraph::UnitigGraph(deBruijnGraph& dbg)
{
	auto&& unbalanced = dbg.getJunctions();
	auto&& out_unbalanced = unbalanced.first;
	auto&& in_unbalanced = unbalanced.second;
	unsigned int index = 0;
	for (auto& v : out_unbalanced)
	{
		std::string curr = v.get_kmer();
		Vertex* source = dbg.getVertex(curr);
		connectUnbalanced(source, &index, curr, dbg, true);
	}
	for (auto& w: in_unbalanced)
	{
		std::string curr = w.get_kmer();
		Vertex* sink = dbg.getVertex(curr);
		connectUnbalanced(sink, &index, curr, dbg, false);
	}
	// DEBUG
	auto&& numV = boost::num_vertices(g_);
	auto&& numE = boost::num_edges(g_);
	std::cerr << "Unitig graph has " << numV << " vertices and " << numE << " edges" << std::endl;

	typedef std::map<UVertex, int> IndexMap;
	IndexMap mapIndex;
	boost::associative_property_map<IndexMap> propmapIndex(mapIndex);
	uvertex_iter vi, vi_end;
	int i = 0;
	for (boost::tie(vi,vi_end) = boost::vertices(g_); vi != vi_end; ++vi)
	    boost::put(propmapIndex,*vi,i++);
	boost::write_graphviz(std::cout, g_, boost::make_label_writer(boost::get(boost::vertex_name_t(),g_)), boost::make_label_writer(boost::get(boost::edge_name_t(),g_)), boost::default_writer(), propmapIndex);
}

UVertex UnitigGraph::addVertex(unsigned int* index, std::string name)
{
	UVertex uv = boost::add_vertex(g_);
	boost::property_map<UGraph, boost::vertex_name_t>::type n = boost::get(boost::vertex_name_t(), g_);
	boost::put(n, uv, name);
	boost::property_map<UGraph, boost::vertex_index1_t>::type idx = boost::get(boost::vertex_index1_t(), g_);
	unsigned int id = (*index)++;
	boost::put(idx, uv, id);
	auto&& ins = std::make_pair(id,uv);
	graph_.insert(ins);
	return uv;
}

void UnitigGraph::connectUnbalanced(Vertex* source, unsigned int* index, std::string curr, deBruijnGraph& dbg, bool forward)
{
	std::vector<char> neigh;
	if (forward)
		neigh = source->get_successors();
	else
		neigh = source->get_predecessors();
	if (!source->is_visited())
	{
		source->visit();
		UVertex uv = addVertex(index, curr);
		source->set_index(*index);
		//dfs for all neighbours
		for (const auto& n : neigh)
		{
			std::string sequence("");
			std::string next;
			if (forward)
				next = curr.substr(1) + n;
			else
				next = n + curr.substr(0,curr.length() - 1);
			sequence += n;
			auto&& nextV = dbg.getVertex(next);
			// what if nextV = 0
			buildEdge(uv, nextV, next, sequence, index, dbg);
		}
	}
	else
	{
		unsigned int idx = source->get_index();
		UVertex uv = graph_[idx];
		for (const auto& n : neigh)
		{
			std::string next;
			if (forward)
				next = curr.substr(1) + n;
			else
				next = n + curr.substr(0,curr.length() - 1);
			auto&& nextV = dbg.getVertex(next);
			if (nextV->is_visited())
				break; // neighbours have been checked
			else
			{
				std::string sequence("");
				sequence += n;
				buildEdge(uv, nextV, next, sequence, index, dbg);
			}
		}
	}
}

UEdge UnitigGraph::buildEdge(UVertex src, Vertex* nextV, std::string next, std::string& sequence, unsigned int* index, deBruijnGraph& dbg)
{
	auto&& succ = nextV->get_successors();
	while (!nextV->is_visited() and succ.size() == 1 and nextV->get_predecessors().size() == 1)
	{
		nextV->visit();
		auto&& c = succ[0];
		Sequence tmp = dbg.getSequence(next);
		if (tmp == next)
			next = next.substr(1) + c;
		else
			next = deBruijnGraph::complement(c) + next.substr(0,next.length() - 1);
		nextV = dbg.getVertex(next);
		succ = nextV->get_successors();
		sequence += c; 
	} //after this we found a junction (visited vertex can only be visited twice if junction)
	if (!nextV->is_visited())
	{
		nextV->visit();
		addVertex(index, next);
		nextV->set_index(*index);
	}
	UEdge e = (boost::add_edge(src, graph_[nextV->get_index()],g_)).first;
	boost::property_map<UGraph, boost::edge_name_t>::type name = boost::get(boost::edge_name_t(), g_);
	boost::put(name,e,sequence);
	return e;
}
