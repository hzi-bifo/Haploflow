#include "UnitigGraph.h"

// so we can use range-based for loops on BGL graphs
namespace std
{
	template <typename A> A begin(const pair<A, A>& s)
	{
		return s.first;
	}

	template <typename A> A end(const pair<A, A>& s)
	{
		return s.second;
	}
}

// constructor of the so-called UnitigGraph
// unifies all simple paths in the deBruijnGraph to a single source->sink path
// all remaining nodes have either indegree != outdegree or indegree == outdegree > 1
UnitigGraph::UnitigGraph(deBruijnGraph& dbg)
{
	// first retrieve all unbalanced vertices in debruijn graph
	auto&& unbalanced = dbg.getJunctions();
	auto&& out_unbalanced = unbalanced.first;
	auto&& in_unbalanced = unbalanced.second;
	unsigned int index = 1;
	// the functions for connecting them does not differ
	for (auto& v : out_unbalanced)
	{
		std::string curr = v.get_kmer();
		Vertex* source = dbg.getVertex(curr);
		index = connectUnbalanced(source, index, curr, dbg);
	}
	for (auto& w: in_unbalanced)
	{
		std::string curr = w.get_kmer();
		Vertex* sink = dbg.getVertex(curr);
		index = connectUnbalanced(sink, index, curr, dbg);
	}
	// DEBUG
	auto&& numV = boost::num_vertices(g_);
	auto&& numE = boost::num_edges(g_);
	std::cerr << "Unitig graph has " << numV << " vertices and " << numE << " edges, starting cleaning" << std::endl;
	
	cleanGraph();
	
	numV = boost::num_vertices(g_);
	numE = boost::num_edges(g_);
	std::cerr << "Finished cleaning" << std::endl;
	std::cerr << "Unitig graph has " << numV << " vertices and " << numE << " edges" << std::endl;

	typedef std::map<UVertex, int> IndexMap;
	IndexMap mapIndex;
	boost::associative_property_map<IndexMap> propmapIndex(mapIndex);
	uvertex_iter vi, vi_end;
	int i = 1;
	int simpletons = 0;
	const auto& cap = boost::get(boost::edge_capacity_t(),g_);
	const auto& name = boost::get(boost::edge_name_t(),g_);
	const auto& vname = boost::get(boost::vertex_name_t(),g_);
	for (boost::tie(vi,vi_end) = boost::vertices(g_); vi != vi_end; ++vi)
	{
		auto indegree = boost::in_degree(*vi,g_);
		auto outdegree = boost::out_degree(*vi,g_);
		if (indegree == 0)
		{
			if (outdegree == 1)
				std::cerr << "Source: " << boost::get(vname,*vi) << " (Vertex " << i << ")" << std::endl;
		}
		else if (outdegree == 0)
		{
			if (indegree == 1)
				std::cerr << "Sink: " << boost::get(vname,*vi) << " (Vertex " << i << ")" <<std::endl;
		}
		else if (indegree == 1 and outdegree == 1)
			simpletons++;
		else
		{
			if (outdegree > 1)
			{
				for (const auto& out : boost::out_edges(*vi,g_))
				{
					auto e = boost::edge(*vi,boost::target(out,g_),g_);
					//std::cerr << boost::get(name,e.first) << ": " << boost::get(cap,e.first) << std::endl;
				}
				//std::cerr << std::endl;
			}
		}
		boost::put(propmapIndex,*vi,i++);
	}
	boost::write_graphviz(std::cout, g_, boost::make_label_writer(boost::get(boost::vertex_name_t(),g_)), boost::make_label_writer(boost::get(boost::edge_name_t(),g_)), boost::default_writer(), propmapIndex);
	//boost::write_graphviz(std::cout, g_, boost::default_writer(), boost::make_label_writer(boost::get(boost::edge_capacity_t(),g_)), boost::default_writer(), propmapIndex);
}

// adds a vertex to the unitig graph: adds it to the boost graph, as well as to the mapping from index to vertex
UVertex UnitigGraph::addVertex(unsigned int index, std::string name)
{
	UVertex uv = boost::add_vertex(g_);
	boost::property_map<UGraph, boost::vertex_name_t>::type n = boost::get(boost::vertex_name_t(), g_);
	boost::put(n, uv, name);
	boost::property_map<UGraph, boost::vertex_index1_t>::type idx = boost::get(boost::vertex_index1_t(), g_);
	boost::put(idx, uv, index);
	auto&& ins = std::make_pair(index,uv);
	graph_.insert(ins);
	return uv;
}

// function for connecting a given source/sink vertex to all its unbalanced successors/predecessors
unsigned int UnitigGraph::connectUnbalanced(Vertex* source, unsigned int index, std::string curr, deBruijnGraph& dbg)
{
	std::vector<char> succ = source->get_successors();
	std::vector<char> pred = source->get_predecessors();
	UVertex uv;
	// if the source is not yet a vertex of the graph, add it
	if (!source->is_visited())
	{
		uv = addVertex(++index, curr);
		source->set_index(index);
		//dfs for all neighbours
	}
	else
	{
		unsigned int idx = source->get_index();
		uv = graph_[idx];
	}
	source->visit(); // make sure the next time we find it we dont add it another time
	index = addNeighbours(curr, succ, pred, dbg, index, uv); // finding the next unbalanced vertices

	return index;
}

// iterating over all neighbours of current node, build the different sequences to the next unbalanced node
unsigned int UnitigGraph::addNeighbours(std::string& curr, const std::vector<char>& succ, const std::vector<char>& pred, deBruijnGraph& dbg, unsigned int index, UVertex& uv)
{
	//bool rc = false;
	// since the source or the sink node has been created before, every edge can at least add one new vertex
	// this means the index can be increased by at most one in every step of the for loop
	// first, find the successing unbalanced vertices
	auto&& currV = dbg.getVertex(curr);
	// check for "real" source/sink property here
	for (const auto& n : succ)
	{
		std::string sequence("");
		std::string next = curr.substr(1) + n;
		Sequence s = dbg.getSequence(next);
		auto&& nextV = dbg.getVertex(next);
		unsigned int coverage = currV->get_out_coverage(n);
		if (s == next)
		{
			sequence += n;
			if (!buildEdge(uv, nextV, next, sequence, index, coverage, dbg))
				index++; // make sure to keep the index intact
		}
		else
		{
			sequence += deBruijnGraph::complement(n);
			if (!buildEdgeReverse(uv, nextV, next, sequence, index, coverage, dbg))
				index++;
		}
	}
	// finding the predecessing unbalanced vertices
	for (const auto& n : pred)
	{
		std::string sequence("");
		std::string prev = n + curr.substr(0,curr.length() - 1);
		Sequence s = dbg.getSequence(prev);
		auto&& nextV = dbg.getVertex(prev);
		unsigned int coverage = currV->get_in_coverage(n);
		if (s == prev)
		{
			sequence += deBruijnGraph::complement(n);
			if (!buildEdgeReverse(uv, nextV, prev, sequence, index, coverage, dbg))
				index++;
		}
		else
		{
			sequence += n;
			if (!buildEdge(uv, nextV, prev, sequence, index, coverage, dbg))
				index++;
		}
	}
	//}
	return index;
}

// go back through the graph until the next unbalanced node is found and add an ("reversed") edge
bool UnitigGraph::buildEdgeReverse(UVertex trg, Vertex* nextV, std::string prev, std::string& sequence, unsigned int index, unsigned int coverage, deBruijnGraph& dbg)
{
	bool visited = true;
	auto&& succ = nextV->get_successors();
	auto&& pred = nextV->get_predecessors();
	// DEBUG
	unsigned int min = 1000000; // TODO
	unsigned int max = 0;
	float avg = coverage;
	unsigned int length = 1;
	// loop until the next unbalanced node is found, which is either visited (has been added) or will be added
	while (!nextV->is_visited() and succ.size() == 1 and pred.size() == 1)
	{
		nextV->visit();
		char c = pred[0];
		unsigned int cov = nextV->get_in_coverage(c); // should be > 0 (assert)
		if (cov < min)
			min = cov;
		if (cov > max)
			max = cov;
		avg += cov;
		length++;
		Sequence tmp = dbg.getSequence(prev);
		if (tmp == prev)
			prev = c + prev.substr(0,prev.length() - 1);
		else
		{
			c = succ[0];
			prev = deBruijnGraph::complement(c) + prev.substr(0,prev.length() - 1);
		}
		nextV = dbg.getVertex(prev);
		pred = nextV->get_predecessors();
		succ = nextV->get_successors();
		sequence += deBruijnGraph::complement(c);
	}
	// if the next vertex has been visited it already is part of the unitiggraph, otherwise add it
	if (!nextV->is_visited())
	{
		nextV->visit();
		addVertex(++index, prev);
		nextV->set_index(index);
		visited = false;
	}
	// this shouldnt happen (it does though) TODO
	else if (succ.size() == 1 and pred.size() == 1)
	{
		return true;
	}
	avg /= float(length); // average
	boost::property_map<UGraph, boost::edge_name_t>::type name = boost::get(boost::edge_name_t(), g_);
	boost::property_map<UGraph, boost::edge_capacity_t>::type cap = boost::get(boost::edge_capacity_t(), g_);
	boost::property_map<UGraph, boost::edge_residual_capacity_t>::type len = boost::get(boost::edge_residual_capacity_t(), g_);
	auto src = graph_[nextV->get_index()];
	auto e = boost::edge(src,trg,g_);
	bool toAdd = true;
	if ((sequence.length() <= dbg.getK() and avg < 50) or avg < 30)
		toAdd = false;
	// if edge has been added or the immediate neighbour is an unbalanced vertex, do not add edge TODO
	if (!(sequence.length() == 1 and e.second and boost::get(name,e.first).length() == 1) and toAdd)
	{
		e = boost::add_edge(src,trg,g_);
		std::reverse(sequence.begin(), sequence.end()); // we add the path from the found node to trg
		boost::put(name,e.first,sequence);
		boost::put(cap,e.first,avg);
		boost::put(len,e.first,sequence.length());
	}
	return visited;
}

// same function like the reverse one, but going forward and finding successors
bool UnitigGraph::buildEdge(UVertex src, Vertex* nextV, std::string next, std::string& sequence, unsigned int index, unsigned int coverage, deBruijnGraph& dbg)
{
	bool visited = true;
	// with a little effort this can be moved inside the while loop for efficiency reasons
	auto&& succ = nextV->get_successors();
	auto&& pred = nextV->get_predecessors();
	// DEBUG, coverage information
	unsigned int min = 1000000; // TODO
	unsigned int max = 0;
	float avg = coverage;
	unsigned int length = 1;
	while (!nextV->is_visited() and succ.size() == 1 and pred.size() == 1)
	{
		nextV->visit();
		char c = succ[0];
		unsigned int cov = nextV->get_out_coverage(c); // should be > 0 (assert)
		if (cov < min)
			min = cov;
		if (cov > max)
			max = cov;
		avg += cov;
		length++;
		Sequence tmp = dbg.getSequence(next);
		if (tmp == next)
			next = next.substr(1) + c;
		else
		{
			c = pred[0];
			next = next.substr(1) + deBruijnGraph::complement(c);
		}
		nextV = dbg.getVertex(next);
		pred = nextV->get_predecessors();
		succ = nextV->get_successors();
		sequence += c; 
	}
	/* 
	if nextV is visited then nextV may either be a junction, in which case it should have been
	added as a vertex to the graph and will receive an edge. Or the path we are starting to build has already
	been found having the target of the path to be found as source. This means we can break now.
	If nextV still isn't visited we found a junction which has not been considered before
	*/
	if (!nextV->is_visited())
	{
		nextV->visit();
		addVertex(++index, next);
		nextV->set_index(index);
		visited = false;
	}
	else if (succ.size() == 1 and pred.size() == 1) // check how this happens
	{
		return true; // path has been found, do not add anything
	}
	avg /= float(length);
	boost::property_map<UGraph, boost::edge_name_t>::type name = boost::get(boost::edge_name_t(), g_);
	boost::property_map<UGraph, boost::edge_capacity_t>::type cap = boost::get(boost::edge_capacity_t(), g_);
	boost::property_map<UGraph, boost::edge_residual_capacity_t>::type len = boost::get(boost::edge_residual_capacity_t(), g_);
	auto trg = graph_[nextV->get_index()];
	auto e = boost::edge(src,trg,g_);
	// error correction, TODO make threshold variable
	bool addEdge = true;
	if ((sequence.length() <= dbg.getK() and avg < 50) or avg < 30) // this path is considered illegal
	{
		/*boost::property_map<UGraph, boost::vertex_name_t>::type vn = boost::get(boost::vertex_name_t(), g_);
		auto& sseq = boost::get(vn, uv);
		std::string sstring = sseq.get_kmer();
		Vertex* v = dbg.getVertex(sstring);
		Vertex* w = dbg.getVertex(next);*/
		addEdge = false;
	}
	if (!(sequence.length() == 1 and e.second and boost::get(name,e.first).length() == 1) and addEdge)
	{
		e = boost::add_edge(src,trg,g_);
		boost::put(name,e.first,sequence);
		boost::put(cap,e.first,avg);
		boost::put(len,e.first,sequence.length());
	}
	return visited;
}

void UnitigGraph::cleanGraph()
{
	boost::graph_traits<UGraph>::vertex_iterator vi, vi_end, next;
	boost::tie(vi, vi_end) = boost::vertices(g_);
	//boost::property_map<UGraph, boost::edge_name_t>::type name = boost::get(boost::edge_name_t(), g_);
  	bool merged = true;
	// as long as some pathes are merged, continue looping
	boost::tie(vi, vi_end) = boost::vertices(g_);
	for (next = vi; vi != vi_end; vi = next) {
		++next;
		auto indegree = boost::in_degree(*vi, g_);
		auto outdegree = boost::out_degree(*vi,g_);
		if (outdegree == 0 and indegree == 0)
		{
			boost::remove_vertex(*vi,g_);
		}
	}
}
