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
UnitigGraph::UnitigGraph(deBruijnGraph& dbg) : cc_(1)
{
	std::cerr << "Building unitig graph from deBruijn graph..." << std::endl;
	clock_t t = clock();
	auto&& junc = dbg.getJunctions();
	unsigned int index = 1;
	auto&& out_unbalanced = junc.first;
	auto&& in_unbalanced = junc.second;
	// starting from the sources, we build the unitig graph
	for (auto& v : out_unbalanced)
	{
		std::string curr = v.get_kmer();
		Vertex* source = dbg.getVertex(curr);
		// make a guess whether we are relevant already (in_coverage should always be 0 though)
		if (!source->is_flagged() and (source->get_total_out_coverage() > 30 or source->get_total_in_coverage() > 30))
		{
			connectUnbalanced(source, &index, curr, dbg); 
			cc_++; // new connected component found
		}
	}
	for (auto& v : in_unbalanced)
	{
		std::string curr = v.get_kmer();
		Vertex* source = dbg.getVertex(curr);
		// make a guess whether we are relevant already (in_coverage should always be 0 though)
		if (!source->is_flagged() and (source->get_total_out_coverage() > 30 or source->get_total_in_coverage() > 30))
		{
			connectUnbalanced(source, &index, curr, dbg); 
			cc_++; // new connected component found
		}
	}
	std::cerr << "Unitig graph succesfully build in " << (clock() - t)/1000000. << " seconds." << std::endl;
	// DEBUG
	auto&& numV = boost::num_vertices(g_);
	auto&& numE = boost::num_edges(g_);
	
	std::cerr << "Unitig graph has " << numV << " vertices and " << numE << " edges, starting cleaning" << std::endl;
	t = clock();

	//TODO move this out
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
	//const auto& name = boost::get(boost::edge_name_t(),g_);
	const auto& vname = boost::get(boost::vertex_name_t(),g_);
	for (boost::tie(vi,vi_end) = boost::vertices(g_); vi != vi_end; ++vi)
	{
		unsigned int indegree = boost::in_degree(*vi,g_);
		unsigned int outdegree = boost::out_degree(*vi,g_);
		
		float total_out = 0.;
		float total_in = 0.;
		
		for (auto&& oe : boost::out_edges(*vi,g_))
		{
			auto&& e = boost::edge(*vi,boost::target(oe,g_),g_);
			total_out += boost::get(cap,e.first).first;
		}

		for (auto&& ie : boost::in_edges(*vi,g_))
		{
			auto&& e = boost::edge(boost::source(ie,g_),*vi,g_);
			total_in += boost::get(cap,e.first).first;
		}
		float epsilon = float(total_out)/total_in; // the gain/loss of this vertex
		unsigned int delta = std::abs(total_out - total_in);

		if (false) //DEBUG
		{
			if (outdegree == 0)
			{
				if (indegree == 1)
					std::cerr << "Sink: " << boost::get(vname,*vi) << " (Vertex " << i << ")" <<std::endl;
				else
					std::cerr << "Unreal sink: " << boost::get(vname,*vi) << " (Vertex " << i << ")" << std::endl;
			}
			else if (indegree == 1 and outdegree == 1)
				simpletons++;
			if ((epsilon > 1.2 or epsilon < 0.8) and delta > 10 and indegree > 0 and outdegree > 0)
			{
				std::cerr << boost::get(vname,*vi) << " (" << epsilon << "(" << total_out << ")/" << i << ")" << std::endl;
			}
		}
		boost::put(propmapIndex,*vi,i++);
	}
	//boost::write_graphviz(std::cout, g_, boost::make_label_writer(boost::get(boost::vertex_name_t(),g_)), boost::make_label_writer(boost::get(boost::edge_name_t(),g_)), boost::default_writer(), propmapIndex);
	boost::write_graphviz(std::cout, g_, boost::default_writer(), boost::make_label_writer(boost::get(boost::edge_residual_capacity_t(),g_)), boost::default_writer(), propmapIndex);
	//boost::write_graphviz(std::cout, g_, boost::make_label_writer(boost::get(boost::vertex_index2_t(),g_)), boost::make_label_writer(boost::get(boost::edge_residual_capacity_t(),g_)), boost::default_writer(), propmapIndex);
}

// adds a vertex to the unitig graph: adds it to the boost graph, as well as to the mapping from index to vertex
UVertex UnitigGraph::addVertex(unsigned int* index, std::string name)
{
	UVertex uv = boost::add_vertex(g_);
	auto n = boost::get(boost::vertex_name_t(), g_);
	boost::put(n, uv, name);
	auto idx = boost::get(boost::vertex_index1_t(), g_);
	(*index)++;
	boost::put(idx, uv, *index);
	auto cc = boost::get(boost::vertex_index2_t(), g_);
	boost::put(cc, uv, cc_);
	auto&& ins = std::make_pair(*index,uv);
	graph_.insert(ins);
	return uv;
}

// function for connecting a given source/sink vertex to all its unbalanced successors/predecessors
void UnitigGraph::connectUnbalanced(Vertex* source, unsigned int* index, std::string curr, deBruijnGraph& dbg)
{
	// FIFO queue for the vertices, so we do not find a reverse complement before the original vertex
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
		bool to_search = false;
		// check if there might be valid paths
		for (const char& c : succ)
		{
			if (junction->get_out_coverage(c) > 30)
			{
				to_search = true; break;
			}
		}
		if (!to_search)
		{
			for (const char& c : pred)
			{
				if (junction->get_in_coverage(c) > 30)
				{
					to_search = true; break;
				}
			}
		}
		UVertex uv;
		// a junction is flagged if it has been on the stack before and does not need to be searched again
		if (junction->is_flagged())
		{
			continue;
		}
		// if the source is not yet a vertex of the graph, add it
		if (!junction->is_visited() and to_search) 
		{
			junction->visit();
			uv = addVertex(index, seq);
			junction->set_index(*index);
			//dfs for all neighbours
		}
		else if (!junction->is_visited()) // this vertex is considered to be an erroneous kmer
		{
			junction->visit();
			junction->flag();
			continue;
		}
		else
		{ //this vertex has already been found as the next junction of some other vertex, do not add again
			unsigned int idx = junction->get_index();
			uv = graph_[idx];
		}
		junction->visit(); // make sure the next time we find it we dont add it another time
		auto&& following = addNeighbours(seq, succ, pred, dbg, index, uv); // finding the next unbalanced vertices
		for (auto v : following)
		{	
			// if no sequence is returned, no vertex was added, so we do not need to continue on this vertex
			if (v.second != "")
				todo.push(v);
		}
	} while (!todo.empty());
}

// iterating over all neighbours of current node, build the different sequences to the next unbalanced node
std::vector<std::pair<Vertex*,std::string> > UnitigGraph::addNeighbours(std::string& curr, const std::vector<char>& succ, const std::vector<char>& pred, deBruijnGraph& dbg, unsigned int* index, UVertex& uv)
{
	std::vector<std::pair<Vertex*,std::string> > following;
	Vertex* currV = dbg.getVertex(curr);
	unsigned int total_out = currV->get_total_out_coverage();
	unsigned int total_in = currV->get_total_in_coverage();
	unsigned int threshold = 30; // make variable, based on e.g. median coverage
	if (currV->is_flagged() or total_out < threshold or total_in < threshold)
	{
		return following; // the neighbours for this vertex have been added
	}
	float epsilon = float(total_out)/total_in; // the gain/loss of this vertex
	Sequence src = dbg.getSequence(curr);
	// this is for checking whether the reverse complement is in the debruijn graph
	bool reverse = false;
	if (src != curr)
	{
		reverse = true; 
	}
	for (const auto& n : succ)
	{
		std::string sequence("");
		std::string next;
		float coverage = currV->get_out_coverage(n);
		float pcov = coverage/float(total_out);// TODO /total_out; // percentage of coverage
		if (!reverse)
		{
			next = curr.substr(1) + n;
			Vertex* nextV = dbg.getVertex(next);
			Sequence s = dbg.getSequence(next);
			sequence += n;
			if (!nextV->is_flagged())
				following.push_back(buildEdge(uv, nextV, next, sequence, index, coverage, pcov, dbg));
		}
		// if we are a reverse complement, we actually want to add the path in reverse order
		else
		{
			next = deBruijnGraph::complement(n) + curr.substr(0,curr.length() - 1);
			Vertex* nextV = dbg.getVertex(next);
			Sequence s = dbg.getSequence(next);
			sequence += curr.back(); // the predecessor points to the current vertex with the last char of curr (by definition)
			if (!nextV->is_flagged())
				following.push_back(buildEdgeReverse(uv, nextV, next, sequence, index, coverage, pcov, dbg));
		}
	}
	// finding the predecessing unbalanced vertices
	for (const auto& n : pred)
	{
		std::string sequence("");
		std::string prev; 
		float coverage = currV->get_in_coverage(n);
		float pcov = coverage/float(total_in); 
		if (!reverse)
		{
			prev = n + curr.substr(0,curr.length() - 1);
			Vertex* nextV = dbg.getVertex(prev);
			Sequence s = dbg.getSequence(prev);
			sequence += curr.back(); // the predecessor points to the current vertex with the last char of curr
			if (!nextV->is_flagged())
				following.push_back(buildEdgeReverse(uv, nextV, prev, sequence, index, coverage, pcov, dbg));
		}
		else
		{
			prev = curr.substr(1) + deBruijnGraph::complement(n);
			Vertex* nextV = dbg.getVertex(prev);
			Sequence s = dbg.getSequence(prev);
			sequence += deBruijnGraph::complement(n);
			if (!nextV->is_flagged())
				following.push_back(buildEdge(uv, nextV, prev, sequence, index, coverage, pcov, dbg));
		}
	}
	currV->flag(); // this vertex is done
	return following;
}

// go back through the graph until the next unbalanced node is found and add an ("reversed") edge
std::pair<Vertex*,std::string> UnitigGraph::buildEdgeReverse(UVertex trg, Vertex* nextV, std::string prev, std::string& sequence, unsigned int* index, float coverage, float pcov, deBruijnGraph& dbg)
{
	float frac = pcov; // the amount of flow for certain char divided by total flow towards this vertex
	auto&& succ = nextV->get_successors();
	auto&& pred = nextV->get_predecessors();
	// DEBUG
	float min = coverage; // TODO
	float max = coverage;
	float avg = coverage;
	unsigned int length = 1;
	auto vn = boost::get(boost::vertex_name_t(), g_);
	char lastchar = boost::get(vn,trg).back(); //char with which we are pointing to trg
	// loop until the next unbalanced node is found, which is either visited (has been added) or will be added
	while (!nextV->is_visited() and succ.size() == 1 and pred.size() == 1)
	{
		nextV->visit();
		Sequence tmp = dbg.getSequence(prev); // check for reverse complimentarity
		unsigned int cov;
		char c;
		lastchar = prev.back(); // this will be added to the sequence
		if (tmp != prev) // we are a reverse complement
		{
			/* if Z<-Y, Y on the complementary strand of Z and Y->X with character c:
			then \overline{Y}<-\overline{X} with character \overline{c} */
			c = deBruijnGraph::complement(succ[0]);
			cov = nextV->get_out_coverage(succ[0]); //percentage will always be 100% on single path
		}
		else
		{
			c = pred[0];
			cov = nextV->get_in_coverage(pred[0]); // see above
		}
		prev = c + prev.substr(0, prev.length() - 1);
		// update the coverage of the path
		if (cov < min)
			min = cov;
		if (cov > max)
			max = cov;
		avg += cov;
		length++; // used for the coverage caluclations of the path later on
		nextV = dbg.getVertex(prev);
		pred = nextV->get_predecessors();
		succ = nextV->get_successors();
		sequence += lastchar;
	}
	avg /= float(length); // average coverage over the path
	// if the next vertex has been visited it already is part of the unitiggraph, otherwise add it
	unsigned int threshold = 50; //TODO variable/based on mean/median/... coverage?
	unsigned int long_thresh = 30; // if the path is longer than e.g. read length, the coverage is allowed to be lower for detecting lowly abundant strains
	
	if (!nextV->is_visited() and ((avg >= threshold or sequence.length() > dbg.getK()) and avg >= long_thresh))//and nextV->get_total_out_coverage() > 30 and nextV->get_total_in_coverage() > 30)
	{
		nextV->visit();
		addVertex(index, prev); // the vertex is new and found to be relevant
		nextV->set_index(*index);
	}
	else if (!nextV->is_visited())
	{
		return std::make_pair(nextV,""); // path has too low coverage
	}
	else if (succ.size() == 1 and pred.size() == 1)
	{
		return std::make_pair(nextV,""); // TODO
		// this vertex has been found on another path
	}
	// if the vertex has not been added, dont try to add an edge to it
	if (nextV->get_index() == 0)
	{
		return std::make_pair(nextV,"");
	}

	auto name = boost::get(boost::edge_name_t(), g_);
	auto cap = boost::get(boost::edge_capacity_t(), g_);
	auto rcap = boost::get(boost::edge_residual_capacity_t(),g_);

	UVertex src = graph_[nextV->get_index()];
	auto e = boost::edge(src,trg,g_);
	bool toAdd = true;
	float min_pass = 5; // minimal coverage allowed on edge
	// path will not be added
	// TODO we should not end up here with avg < threshold etc.
	if ((sequence.length() <= dbg.getK() and avg < threshold) or avg < long_thresh)
	{
		toAdd = false;
	}
	else if (min < min_pass)
	{
		// the path has a high enough average, but some edges on it had a low coverage
		// TODO indel detected?
		//toAdd = false; // this is interesting!
	}
	// if edge has been added or the immediate neighbour is an unbalanced vertex, do not add edge TODO
	if ((!e.second or (e.second and (boost::get(name,e.first)) != sequence)) and toAdd)
	{
		e = boost::add_edge(src,trg,g_);
		std::reverse(sequence.begin(), sequence.end()); // we add the path from the found node to trg, the sequence was added in reverse order
		boost::put(name,e.first,sequence);
		boost::put(cap,e.first,std::make_pair(avg,frac));
		Sequence s = dbg.getSequence(prev);
		if (s != prev)
		{
			float in = nextV->get_in_coverage(deBruijnGraph::complement(lastchar));
			float tot_in = nextV->get_total_in_coverage();
			boost::put(rcap,e.first,in/tot_in);
		}
		else
		{
			float out = nextV->get_out_coverage(lastchar);
			float tot_out = nextV->get_total_out_coverage();
			boost::put(rcap,e.first,out/tot_out);
		}
	}
	return std::make_pair(nextV,prev);
}

// same function like the reverse one, but going forward and finding successors
std::pair<Vertex*,std::string> UnitigGraph::buildEdge(UVertex src, Vertex* nextV, std::string next, std::string& sequence, unsigned int* index, float coverage, float pcov, deBruijnGraph& dbg)
{
	float frac = pcov;// the amount of flow for certain char divided by total flow towards this vertex
	// with a little effort this can be moved inside the while loop for efficiency reasons
	auto&& succ = nextV->get_successors();
	auto&& pred = nextV->get_predecessors();
	// DEBUG, coverage information
	float min = coverage; // TODO
	float max = coverage;
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
	unsigned int threshold = 50;
	unsigned int long_thresh = 30; // if the path is longer than e.g. read length, the coverage is allowed to be lower for detecting lowly abundant strains
	if (!nextV->is_visited() and ((avg >= threshold or sequence.length() > dbg.getK()) and avg >= long_thresh)) // TODO dbg.getK() should e.g. be dbg.readLength() 
	{
		nextV->visit();
		addVertex(index, next);
		nextV->set_index(*index);
	}
	else if (!nextV->is_visited())
	{
		return std::make_pair(nextV,"");
	}
	else if (succ.size() == 1 and pred.size() == 1) // check else-if case
	{
		//TODO
		return std::make_pair(nextV,""); // path has been found, do not add anything
	}
	// this vertex has been found and considered irrelevant because of insufficient coverage
	if (nextV->get_index() == 0)
	{
		return std::make_pair(nextV,"");
	}
	
	auto name = boost::get(boost::edge_name_t(), g_);
	auto cap = boost::get(boost::edge_capacity_t(), g_);
	auto rcap = boost::get(boost::edge_residual_capacity_t(),g_);
	
	UVertex trg = graph_[nextV->get_index()];
	auto e = boost::edge(src,trg,g_);
	// error correction, TODO make threshold variable
	bool toAdd = true;
	float min_pass = 5; // minimal coverage allowed on edge
	if ((sequence.length() <= dbg.getK() and avg < threshold) or avg < long_thresh) // this path is considered illegal
	{
		/*boost::property_map<UGraph, boost::vertex_name_t>::type vn = boost::get(boost::vertex_name_t(), g_);
		auto& sseq = boost::get(vn, uv);
		std::string sstring = sseq.get_kmer();
		Vertex* v = dbg.getVertex(sstring);
		Vertex* w = dbg.getVertex(next);*/
		toAdd = false;
	}
	else if (min < min_pass)
	{
		// TODO
	}
	if ((!e.second or (e.second and (boost::get(name,e.first)) != sequence)) and toAdd)
	{
		e = boost::add_edge(src,trg,g_);
		boost::put(name,e.first,sequence);
		boost::put(cap,e.first,std::make_pair(avg,frac));
		boost::put(rcap,e.first,frac);
	}
	return std::make_pair(nextV,next);
}

// the graph might contain some unconnected vertices, clean up
// TODO maybe some more cleaning (e.g. long, simple paths)
void UnitigGraph::cleanGraph()
{
	boost::graph_traits<UGraph>::vertex_iterator vi, vi_end, next;
	//boost::property_map<UGraph, boost::vertex_name_t>::type name = boost::get(boost::vertex_name_t(), g_);
	const auto& cap = boost::get(boost::edge_capacity_t(),g_);
	const auto& rcap = boost::get(boost::edge_residual_capacity_t(),g_);
	const auto& name = boost::get(boost::edge_name_t(),g_);
	boost::tie(vi, vi_end) = boost::vertices(g_);
	for (next = vi; vi != vi_end; vi = next)
	{
		++next;
		unsigned int indegree = boost::in_degree(*vi, g_);
		unsigned int outdegree = boost::out_degree(*vi,g_);
		
		//if (false) // this might be too strict, capacity information is not really preserved
		// if in and outdegree is 1, we are on a simple path and can contract again
		if (outdegree == 1 and indegree == 1)
		{
			auto&& ie = boost::in_edges(*vi,g_);
			auto&& oe = boost::out_edges(*vi,g_);
			auto&& new_source = boost::source(*ie.first,g_);
			auto&& new_target = boost::target(*oe.first,g_);
			auto&& e = boost::edge(new_source,*vi,g_);
			std::string seq = boost::get(name,e.first);
			unsigned int w = seq.length();
			float capacity = boost::get(cap,e.first).first * w;
			e = boost::edge(*vi,new_target,g_);
			seq += boost::get(name,e.first);
			capacity += (boost::get(cap,e.first).first * (seq.length() - w));
			capacity /= seq.length();
			e = boost::add_edge(new_source,new_target,g_);
			boost::put(name,e.first,seq);
			boost::put(cap,e.first,std::make_pair(capacity,capacity)); // TODO
			boost::put(rcap,e.first,capacity);
			boost::clear_vertex(*vi,g_);
			boost::remove_vertex(*vi,g_);
		}
	}
	// DEBUG
	boost::graph_traits<UGraph>::edge_iterator ei, ei_end;
	boost::tie(ei, ei_end) = boost::edges(g_);
	for (; ei != ei_end; ei++)
	{
		float c = boost::get(cap,*ei).first;
		boost::put(rcap,*ei,c);
	}
	boost::tie(vi, vi_end) = boost::vertices(g_);
	for (next = vi; vi != vi_end; vi = next) {
		++next;
		unsigned int indegree = boost::in_degree(*vi, g_);
		unsigned int outdegree = boost::out_degree(*vi,g_);
		if (outdegree == 0 and indegree == 0)
		{
			boost::remove_vertex(*vi,g_);
		}
	}
}

// calculates the flows and corresponding paths through the graph
void UnitigGraph::calculateFlow()
{
	
}
