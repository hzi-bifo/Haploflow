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
	std::cerr << "Building unitig graph from deBruijn graph..." << std::endl;
	clock_t t = clock();
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
		if (!source->is_flagged() and source->get_total_out_coverage() > 30 and source->get_total_in_coverage() > 30)
		{
			connectUnbalanced(source, &index, curr, dbg);
		}
	}
	for (auto& w: in_unbalanced)
	{
		break;
		std::string curr = w.get_kmer();
		Vertex* sink = dbg.getVertex(curr);
		if (!sink->is_flagged() and sink->get_total_out_coverage() > 30 and sink->get_total_in_coverage() > 30)
		{
			connectUnbalanced(sink, &index, curr, dbg);
		}
	}
	std::cerr << "Unitig graph succesfully build in " << (clock() - t)/1000000. << " seconds." << std::endl;
	// DEBUG
	auto&& numV = boost::num_vertices(g_);
	auto&& numE = boost::num_edges(g_);
	
	std::cerr << "Unitig graph has " << numV << " vertices and " << numE << " edges, starting cleaning" << std::endl;
	t = clock();

	cleanGraph();
	
	numV = boost::num_vertices(g_);
	numE = boost::num_edges(g_);
	std::cerr << "Finished cleaning: Unitig graph has " << numV << " vertices and " << numE << " edges" << std::endl;
	std::cerr << "Cleaning took " << (clock() - t)/1000000. << " seconds." << std::endl;

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
		if (false)
		{
			if (indegree == 0)
			{
				if (outdegree == 1)
					std::cerr << "Source: " << boost::get(vname,*vi) << " (Vertex " << i << ")" << std::endl;
				//else
				//	std::cerr << "Unreal source: " << boost::get(vname,*vi) << " (Vertex " << i << ")" << std::endl;
			}
			else if (outdegree == 0)
			{
				if (indegree == 1)
					std::cerr << "Sink: " << boost::get(vname,*vi) << " (Vertex " << i << ")" <<std::endl;
				//else
				//	std::cerr << "Unreal sink: " << boost::get(vname,*vi) << " (Vertex " << i << ")" << std::endl;
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
						std::cerr << boost::get(name,e.first) << ": " << boost::get(cap,e.first) << std::endl;
					}
					std::cerr << "___" << std::endl;
				}
				if (indegree > 1)
				{
					for (const auto& out : boost::in_edges(*vi,g_))
					{
						auto e = boost::edge(boost::source(out,g_),*vi,g_);
						std::cerr << boost::get(name,e.first) << ": " << boost::get(cap,e.first) << std::endl;
					}
					std::cerr << "____" <<std::endl;
				}
			}
		}
		boost::put(propmapIndex,*vi,i++);
	}
	boost::write_graphviz(std::cout, g_, boost::make_label_writer(boost::get(boost::vertex_name_t(),g_)), boost::make_label_writer(boost::get(boost::edge_name_t(),g_)), boost::default_writer(), propmapIndex);
	//boost::write_graphviz(std::cout, g_, boost::default_writer(), boost::make_label_writer(boost::get(boost::edge_capacity_t(),g_)), boost::default_writer(), propmapIndex);
}

// adds a vertex to the unitig graph: adds it to the boost graph, as well as to the mapping from index to vertex
UVertex UnitigGraph::addVertex(unsigned int* index, std::string name)
{
	UVertex uv = boost::add_vertex(g_);
	boost::property_map<UGraph, boost::vertex_name_t>::type n = boost::get(boost::vertex_name_t(), g_);
	boost::put(n, uv, name);
	boost::property_map<UGraph, boost::vertex_index1_t>::type idx = boost::get(boost::vertex_index1_t(), g_);
	(*index)++;
	boost::put(idx, uv, *index);
	auto&& ins = std::make_pair(*index,uv);
	graph_.insert(ins);
	return uv;
}

// function for connecting a given source/sink vertex to all its unbalanced successors/predecessors
void UnitigGraph::connectUnbalanced(Vertex* source, unsigned int* index, std::string curr, deBruijnGraph& dbg)
{
	std::queue<std::pair<Vertex*,std::string> > todo;
	todo.push(std::make_pair(source,curr));
	do 
	{
		auto next = todo.front();
		todo.pop();
		Vertex* junction = next.first;
		std::string seq = next.second;
		std::vector<char> succ = junction->get_successors();
		std::vector<char> pred = junction->get_predecessors();
		UVertex uv;
		// if the source is not yet a vertex of the graph, add it
		if (junction->is_flagged())
		{
			continue;
		}
		if (!junction->is_visited() and junction->get_total_out_coverage() > 30 and junction->get_total_in_coverage() > 30)
		{
			uv = addVertex(index, curr);
			junction->set_index(*index);
			//dfs for all neighbours
		}
		else if (!junction->is_visited())
		{
			junction->visit();
			junction->flag();
			continue;
		}
		else
		{
			unsigned int idx = junction->get_index();
			uv = graph_[idx];
		}
		junction->visit(); // make sure the next time we find it we dont add it another time
		auto&& following = addNeighbours(seq, succ, pred, dbg, index, uv); // finding the next unbalanced vertices
		for (auto v : following)
		{
			todo.push(v);
		}
	} while (!todo.empty());
}

// iterating over all neighbours of current node, build the different sequences to the next unbalanced node
std::vector<std::pair<Vertex*,std::string> > UnitigGraph::addNeighbours(std::string& curr, const std::vector<char>& succ, const std::vector<char>& pred, deBruijnGraph& dbg, unsigned int* index, UVertex& uv)
{
	//bool rc = false;
	// since the source or the sink node has been created before, every edge can at least add one new vertex
	// this means the index can be increased by at most one in every step of the for loop
	// first, find the successing unbalanced vertices
	std::vector<std::pair<Vertex*,std::string> > following;
	auto&& currV = dbg.getVertex(curr);
	if (currV->is_flagged() or currV->get_total_out_coverage() < 30 or currV->get_total_in_coverage() < 30)
	{
		return following; // the neighbours for this vertex have been added
	}
	Sequence src = dbg.getSequence(curr);
	bool reverse = false;
	if (src != curr)
	{
		reverse = true; 
	}
	// check for "real" source/sink property here?
	for (const auto& n : succ)
	{
		std::string sequence("");
		std::string next;
		unsigned int coverage = currV->get_out_coverage(n);
		if (!reverse)
		{
			next = curr.substr(1) + n;
			auto&& nextV = dbg.getVertex(next);
			sequence += n;
			if (!nextV->is_flagged())
				following.push_back(buildEdge(uv, nextV, next, sequence, index, coverage, dbg));
		}
		else
		{
			next = deBruijnGraph::complement(n) + curr.substr(0,curr.length() - 1);
			auto&& nextV = dbg.getVertex(next);
			sequence += curr.back();//deBruijnGraph::complement(n);
			if (!nextV->is_flagged())
				following.push_back(buildEdgeReverse(uv, nextV, next, sequence, index, coverage, dbg));
		}
	}
	// finding the predecessing unbalanced vertices
	for (const auto& n : pred)
	{
		std::string sequence("");
		std::string prev; 
		unsigned int coverage = currV->get_in_coverage(n);
		if (!reverse)
		{
			prev = n + curr.substr(0,curr.length() - 1);
			auto&& nextV = dbg.getVertex(prev);
			sequence += curr.back();//n;
			if (!nextV->is_flagged())
				following.push_back(buildEdgeReverse(uv, nextV, prev, sequence, index, coverage, dbg));
		}
		else
		{
			prev = curr.substr(1) + deBruijnGraph::complement(n);
			auto&& nextV = dbg.getVertex(prev);
			sequence += deBruijnGraph::complement(n);
			if (!nextV->is_flagged())
				following.push_back(buildEdge(uv, nextV, prev, sequence, index, coverage, dbg));
		}
	}
	currV->flag(); // this vertex is done
	return following;
}

// go back through the graph until the next unbalanced node is found and add an ("reversed") edge
std::pair<Vertex*,std::string> UnitigGraph::buildEdgeReverse(UVertex trg, Vertex* nextV, std::string prev, std::string& sequence, unsigned int* index, unsigned int coverage, deBruijnGraph& dbg)
{
	auto&& succ = nextV->get_successors();
	auto&& pred = nextV->get_predecessors();
	// DEBUG
	unsigned int min = coverage; // TODO
	unsigned int max = coverage;
	float avg = coverage;
	unsigned int length = 1;
	// loop until the next unbalanced node is found, which is either visited (has been added) or will be added
	while (!nextV->is_visited() and succ.size() == 1 and pred.size() == 1)
	{
		nextV->visit();
		Sequence tmp = dbg.getSequence(prev); // check for reverse complimentarity
		char c;
		unsigned int cov;
		char lastchar = prev.back();
		if (tmp != prev) // we are a reverse complement
		{
			/* if Z<-Y, Y on the complementary strand of Z and Y->X with character c:
			then \overline{Y}<-\overline{X} with character \overline{c} */
			c = deBruijnGraph::complement(succ[0]);
			cov = nextV->get_out_coverage(succ[0]);
		}
		else
		{
			c = pred[0];
			cov = nextV->get_in_coverage(pred[0]);
		}
		prev = c + prev.substr(0, prev.length() - 1);
		if (cov < min)
			min = cov;
		if (cov > max)
			max = cov;
		avg += cov;
		length++;
		nextV = dbg.getVertex(prev);
		pred = nextV->get_predecessors();
		succ = nextV->get_successors();
		sequence += lastchar;
	}
	avg /= float(length); // average
	// if the next vertex has been visited it already is part of the unitiggraph, otherwise add it
	if (!nextV->is_visited() and nextV->get_total_out_coverage() > 30 and nextV->get_total_in_coverage() > 30)
	{
		nextV->visit();
		addVertex(index, prev);
		nextV->set_index(*index);
	}
	else if (!nextV->is_visited())
	{
		nextV->visit();
		nextV->flag();
		return std::make_pair(nextV,prev);
	}
	// this shouldnt happen (it does though) TODO
	else if (succ.size() == 1 and pred.size() == 1)
	{
		return std::make_pair(nextV,prev);
	}
	if (nextV->get_index() == 0)
		return std::make_pair(nextV,prev);
	boost::property_map<UGraph, boost::edge_name_t>::type name = boost::get(boost::edge_name_t(), g_);
	boost::property_map<UGraph, boost::edge_capacity_t>::type cap = boost::get(boost::edge_capacity_t(), g_);
	boost::property_map<UGraph, boost::edge_residual_capacity_t>::type len = boost::get(boost::edge_residual_capacity_t(), g_);
	auto src = graph_[nextV->get_index()];
	auto e = boost::edge(src,trg,g_);
	bool toAdd = true;
	if ((sequence.length() <= dbg.getK() and avg < 50) or avg < 30)
	{
		toAdd = false;
	}
	else if (min < 5)
	{
		//boost::property_map<UGraph, boost::vertex_name_t>::type vn = boost::get(boost::vertex_name_t(), g_);
		//std::cerr << boost::get(vn,src) << " - " << boost::get(vn,trg) << " (" << sequence << ")" << std::endl;
		//toAdd = false; // this is interesting!
	}
	// if edge has been added or the immediate neighbour is an unbalanced vertex, do not add edge TODO
	if ((!e.second or (e.second and (boost::get(name,e.first)) != sequence)) and toAdd)
	{
		e = boost::add_edge(src,trg,g_);
		std::reverse(sequence.begin(), sequence.end()); // we add the path from the found node to trg
		boost::put(name,e.first,sequence);
		boost::put(cap,e.first,avg);
		boost::put(len,e.first,sequence.length());
	}
	return std::make_pair(nextV,prev);
}

// same function like the reverse one, but going forward and finding successors
std::pair<Vertex*,std::string> UnitigGraph::buildEdge(UVertex src, Vertex* nextV, std::string next, std::string& sequence, unsigned int* index, unsigned int coverage, deBruijnGraph& dbg)
{
	// with a little effort this can be moved inside the while loop for efficiency reasons
	auto&& succ = nextV->get_successors();
	auto&& pred = nextV->get_predecessors();
	// DEBUG, coverage information
	unsigned int min = coverage; // TODO
	unsigned int max = coverage;
	float avg = coverage;
	unsigned int length = 1;
	while (!nextV->is_visited() and succ.size() == 1 and pred.size() == 1)
	{
		nextV->visit();
		Sequence tmp = dbg.getSequence(next);
		unsigned int cov;
		char c;
		if (tmp != next) // we are a reverse complement
		{
			/* if X->Y, Y on the complementary strand of X and Y<-Z with character c:
			then \overline{Y}->\overline{Z} with character \overline{c} */
			c = deBruijnGraph::complement(pred[0]);
			cov = nextV->get_in_coverage(pred[0]);
		}
		else
		{
			c = succ[0];
			cov = nextV->get_out_coverage(succ[0]);
		}
		if (cov < min)
			min = cov;
		if (cov > max)
			max = cov;
		avg += cov;
		length++;
		next = next.substr(1) + c;
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
	avg /= float(length);
	if (!nextV->is_visited() and nextV->get_total_out_coverage() > 30 and nextV->get_total_in_coverage() > 30) 
	{
		nextV->visit();
		addVertex(index, next);
		nextV->set_index(*index);
	}
	else if (!nextV->is_visited())
	{
		nextV->visit();
		nextV->flag();
		return std::make_pair(nextV,next);
	}
	else if (succ.size() == 1 and pred.size() == 1) // check else-if case
	{
		//TODO
		return std::make_pair(nextV,next); // path has been found, do not add anything
	}
	// this vertex has been found and considered irrelevant because of insufficient coverage
	if (nextV->get_index() == 0)
		return std::make_pair(nextV,next);
	boost::property_map<UGraph, boost::edge_name_t>::type name = boost::get(boost::edge_name_t(), g_);
	boost::property_map<UGraph, boost::edge_capacity_t>::type cap = boost::get(boost::edge_capacity_t(), g_);
	boost::property_map<UGraph, boost::edge_residual_capacity_t>::type len = boost::get(boost::edge_residual_capacity_t(), g_);
	auto trg = graph_[nextV->get_index()];
	auto e = boost::edge(src,trg,g_);
	// error correction, TODO make threshold variable
	bool toAdd = true;
	if ((sequence.length() <= dbg.getK() and avg < 50) or avg < 30) // this path is considered illegal
	{
		/*boost::property_map<UGraph, boost::vertex_name_t>::type vn = boost::get(boost::vertex_name_t(), g_);
		auto& sseq = boost::get(vn, uv);
		std::string sstring = sseq.get_kmer();
		Vertex* v = dbg.getVertex(sstring);
		Vertex* w = dbg.getVertex(next);*/
		toAdd = false;
	}
	else if (min < 5)
	{
		//boost::property_map<UGraph, boost::vertex_name_t>::type vn = boost::get(boost::vertex_name_t(), g_);
		//std::cerr << boost::get(vn,src) << " - " << boost::get(vn,trg) << " (" << sequence << ")" << std::endl;
	}
	if ((!e.second or (e.second and (boost::get(name,e.first)) != sequence)) and toAdd)
	{
		e = boost::add_edge(src,trg,g_);
		boost::put(name,e.first,sequence);
		boost::put(cap,e.first,avg);
		boost::put(len,e.first,sequence.length());
	}
	return std::make_pair(nextV,next);
}

void UnitigGraph::cleanGraph()
{
	boost::graph_traits<UGraph>::vertex_iterator vi, vi_end, next;
	boost::tie(vi, vi_end) = boost::vertices(g_);
	//boost::property_map<UGraph, boost::vertex_name_t>::type name = boost::get(boost::vertex_name_t(), g_);
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

void UnitigGraph::calculateFlow() const
{

}
