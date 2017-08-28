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

// TODO this is only temporary!
#define THRESHOLD 0
#define LONG_THRESH 0 // threshold if path is not long enough
#define CONTIG_THRESH 0 // if contigs are longer than this they are produced
#define FRAC_THRESH 0.08 // threshold of total coverage from which something is treated significant

// constructor of the so-called UnitigGraph
// unifies all simple paths in the deBruijnGraph to a single source->sink path
// all remaining nodes have either indegree != outdegree or indegree == outdegree > 1
UnitigGraph::UnitigGraph(deBruijnGraph& dbg) : cc_(1)
{
	std::cerr << "deBruijnGraph has " << dbg.getSize() << " vertices" << std::endl;
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
		// make a guess whether we are relevant already 
		if (!source->is_flagged() and (source->get_total_out_coverage() > THRESHOLD or source->get_total_in_coverage() > THRESHOLD))
		{
			connectUnbalanced(source, &index, curr, dbg); 
			cc_++; // new connected component found
		}
	}
	for (auto& v : in_unbalanced)
	{
		std::string curr = v.get_kmer();
		Vertex* source = dbg.getVertex(curr);
		if (!source->is_flagged() and (source->get_total_out_coverage() > THRESHOLD or source->get_total_in_coverage() > THRESHOLD))
		{
			connectUnbalanced(source, &index, curr, dbg); 
			cc_++; // new connected component found
		}
	}
	std::cerr << "Unitig graph succesfully build in " << (clock() - t)/1000000. << " seconds." << std::endl;
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
	boost::put(cc, uv, 0); // sccs are set later
	//boost::put(cc, uv, cc_);
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
	do //we don't necessarily need a do-while anymore
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
			if (junction->get_out_coverage(c) > THRESHOLD)
			{
				to_search = true; break;
			}
		}
		if (!to_search)
		{
			for (const char& c : pred)
			{
				if (junction->get_in_coverage(c) > THRESHOLD)
				{
					to_search = true; break;
				}
			}
		}
		UVertex uv;
		// a junction is flagged if it has been in the queue before and does not need to be searched again
		if (junction->is_flagged())
		{
			continue;
		}
		// if the source is not yet a vertex of the graph, add it
		if (!junction->is_visited() and to_search) 
		{
			junction->visit();
			uv = addVertex(index, seq);
			junction->index = *index;
		}
		else if (!junction->is_visited()) // this vertex is considered to be an erroneous kmer
		{
			junction->visit();
			junction->flag();
			continue;
		}
		else
		{ //this vertex has already been found as the next junction of some other vertex, do not add again
			unsigned int idx = junction->index;
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
	if (currV->is_flagged() or total_out < THRESHOLD or total_in < THRESHOLD)
	{
		return following; // the neighbours for this vertex have been added
	}
	float epsilon = float(total_out)/total_in; // the gain/loss of this vertex
	Sequence src = *dbg.getSequence(curr);
	// this is for checking whether the reverse complement is in the debruijn graph
	bool reverse = (src != curr); // true, if reverse complement, i.e. src and curr are not the same
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
			Sequence s = *dbg.getSequence(next);
			sequence += n;
			if (!nextV->is_flagged())
				following.push_back(buildEdge(uv, nextV, next, sequence, index, coverage, pcov, dbg));
		}
		// if we are a reverse complement, we actually want to add the path in reverse order
		else
		{
			next = deBruijnGraph::complement(n) + curr.substr(0,curr.length() - 1);
			Vertex* nextV = dbg.getVertex(next);
			Sequence s = *dbg.getSequence(next);
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
			Sequence s = *dbg.getSequence(prev);
			sequence += curr.back(); // the predecessor points to the current vertex with the last char of curr
			if (!nextV->is_flagged())
				following.push_back(buildEdgeReverse(uv, nextV, prev, sequence, index, coverage, pcov, dbg));
		}
		else
		{
			prev = curr.substr(1) + deBruijnGraph::complement(n);
			Vertex* nextV = dbg.getVertex(prev);
			Sequence s = *dbg.getSequence(prev);
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
		Sequence tmp = *dbg.getSequence(prev); // check for reverse complimentarity
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
	
	if (!nextV->is_visited() and ((avg >= LONG_THRESH or sequence.length() > dbg.getK()) and avg >= THRESHOLD))//and nextV->get_total_out_coverage() > 30 and nextV->get_total_in_coverage() > 30)
	{
		nextV->visit();
		addVertex(index, prev); // the vertex is new and found to be relevant
		nextV->index = *index;
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
	if (nextV->index == 0)
	{
		return std::make_pair(nextV,"");
	}

	auto name = boost::get(boost::edge_name_t(), g_);
	auto cap = boost::get(boost::edge_capacity_t(), g_);
	auto visit = boost::get(boost::edge_index_t(), g_);

	UVertex src = graph_[nextV->index];
	auto e = boost::edge(src,trg,g_);
	bool toAdd = true;
	float min_pass = 5; // minimal coverage allowed on edge
	// path will not be added
	if ((sequence.length() <= dbg.getK() and avg < LONG_THRESH) or avg < THRESHOLD)
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
		std::string old_name = boost::get(name, e.first);
		if (e.second) // TODO v-S->w-T->v is treated like v<-S-w<-T-v (should be ST self-loop and TS self-loop)
		{
			boost::put(name,e.first,old_name + sequence);
		}
		else
		{
			boost::put(name,e.first,sequence);
		}
		boost::put(cap,e.first,avg);
		boost::put(visit,e.first,false);
	}
	return std::make_pair(nextV,prev);
}

// same function like the reverse one, but going forward and finding successors
std::pair<Vertex*,std::string> UnitigGraph::buildEdge(UVertex src, Vertex* nextV, std::string next, std::string& sequence, unsigned int* index, float coverage, float pcov, deBruijnGraph& dbg)
{
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
		Sequence tmp = *dbg.getSequence(next);
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
	if (!nextV->is_visited() and ((avg >= LONG_THRESH or sequence.length() > dbg.getK()) and avg >= THRESHOLD)) // TODO dbg.getK() should e.g. be dbg.readLength() 
	{
		nextV->visit();
		addVertex(index, next);
		nextV->index = *index;
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
	if (nextV->index == 0)
	{
		return std::make_pair(nextV,"");
	}
	
	auto name = boost::get(boost::edge_name_t(), g_);
	auto cap = boost::get(boost::edge_capacity_t(), g_);
	auto visit = boost::get(boost::edge_index_t(), g_);
	
	UVertex trg = graph_[nextV->index];
	auto e = boost::edge(src,trg,g_);
	bool toAdd = true;
	float min_pass = 5; // minimal coverage allowed on edge
	if ((sequence.length() <= dbg.getK() and avg < LONG_THRESH) or avg < THRESHOLD) // this path is considered illegal
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
		auto old_name = boost::get(name, e.first);
		if (e.second) // TODO v-S->w-T->v is treated like v<-S-w<-T-v (should be ST self-loop and TS self-loop)
		{
			boost::put(name,e.first,sequence + old_name);
		}
		else
		{
			boost::put(name,e.first,sequence);
		}
		boost::put(cap,e.first,avg);
		boost::put(visit,e.first,false);
	}
	return std::make_pair(nextV,next);
}

// contracts all simple paths in graph to a single source-sink connection
void UnitigGraph::contractPaths()
{
	auto cap = boost::get(boost::edge_capacity_t(),g_);
	auto name = boost::get(boost::edge_name_t(),g_);
	auto visit = boost::get(boost::edge_index_t(), g_);
	
	boost::graph_traits<UGraph>::vertex_iterator vi, vi_end, next;
	boost::tie(vi, vi_end) = boost::vertices(g_);
	for (next = vi; vi != vi_end; vi = next)
	{
		++next;
		unsigned int indegree = boost::in_degree(*vi, g_);
		unsigned int outdegree = boost::out_degree(*vi,g_);
		
		// if in and outdegree is 1, we are on a simple path and can contract again TODO coverage calculation
		if (outdegree == 1 and indegree == 1)
		{
			auto&& ie = boost::in_edges(*vi,g_);
			auto&& oe = boost::out_edges(*vi,g_);
			auto&& new_source = boost::source(*ie.first,g_);
			auto&& new_target = boost::target(*oe.first,g_);
			auto&& e = boost::edge(new_source,*vi,g_);
			std::string seq = boost::get(name,e.first);
			unsigned int w = seq.length();
			float capacity = boost::get(cap,e.first) * w;
			e = boost::edge(*vi,new_target,g_);
			seq += boost::get(name,e.first); // append the sequence
			capacity += (boost::get(cap,e.first) * (seq.length() - w));
			capacity /= seq.length(); // currently using the average coverage on the contracted path
			e = boost::add_edge(new_source,new_target,g_);
			boost::put(name,e.first,seq);
			boost::put(cap,e.first,capacity);
			boost::put(visit,e.first,false);
			boost::clear_vertex(*vi,g_);
			boost::remove_vertex(*vi,g_);
		}
	}
}

// checks for a path source/sink to junction whether it is assumed to be siginifcant and returns the path if it isnt
std::vector<UVertex> UnitigGraph::flagDanglingEnd(UVertex& v, bool source)
{
	auto cap = boost::get(boost::edge_capacity_t(),g_);
	
	std::vector<UVertex> pathToJunction{v};
	auto vertex = v;

	int indegree = boost::in_degree(v,g_);
	int outdegree = boost::out_degree(v,g_);
	int coverage = 
		source ? boost::get(cap, *(boost::out_edges(vertex, g_).first)) // coverage of outgoing edge
				: boost::get(cap, *(boost::in_edges(vertex, g_).first)); //coverage of ingoing edge

	while (outdegree <= 1 and indegree <= 1)
	{
		if ((source and outdegree == 0) or (!source and indegree == 0)) // simple path (maybe they will be deleted later on)
		{
			std::vector<UVertex> noDelete{};
			return noDelete;
		}
		auto e = source ? *(boost::out_edges(vertex, g_).first) : *(boost::in_edges(vertex, g_).first);
		vertex = source ? boost::target(e,g_) : boost::source(e,g_);
		pathToJunction.push_back(vertex);
		outdegree = boost::out_degree(vertex,g_);
		indegree = boost::in_degree(vertex,g_);
	} // follow path until next junction
	pathToJunction.pop_back(); // remove the junction itself again, so it doesnt get deleted
	if ((source and indegree > 1 and outdegree == 1) or (!source and indegree == 1 and outdegree > 1)) // next junction
	{
		if (pathToJunction.size() < LONG_THRESH) // path is long enough (TODO make this parameter as well?)
		{
			return pathToJunction;
		}
		float total_coverage = 0.;
		if (source)
		{
			for (const auto& edge : boost::in_edges(vertex, g_))
			{
				total_coverage += boost::get(cap, edge);
			}
		}
		else
		{
			for (const auto& edge : boost::out_edges(vertex, g_))
			{
				total_coverage += boost::get(cap, edge);
			}
		}
		if (coverage/total_coverage < std::min(4*FRAC_THRESH, 0.5)) //TODO since the path is really short, we can be quite strict!
		{
			return pathToJunction;
		}
		else // fraction is okay, we are considered a real source/sink
		{
			std::vector<UVertex> noDelete{};
			return noDelete;
		}
	}
	else
	{
		std::vector<UVertex> noDelete{};
		return noDelete;
	}
}

void UnitigGraph::filterTerminals()
{
	auto vertices = boost::vertices(g_);
	std::vector<UVertex> toDelete;

	for (auto& v : vertices)
	{
		int indegree = boost::in_degree(v,g_);
		int outdegree = boost::out_degree(v,g_);
		if (indegree == 0 or outdegree == 0)
		{
			if (indegree == 0 and outdegree <= 1) // TODO if we are a real junction ourselves, we are assumed to be valid
			{
				std::vector<UVertex> pathToJunction = flagDanglingEnd(v, true);
				toDelete.insert(toDelete.end(), pathToJunction.begin(), pathToJunction.end());
			}
			else if (indegree <= 1 and outdegree == 0)
			{
				std::vector<UVertex> pathToJunction = flagDanglingEnd(v, false);
				toDelete.insert(toDelete.end(), pathToJunction.begin(), pathToJunction.end());
			}
		}
	}
	for (auto& v : toDelete)
	{
		boost::clear_vertex(v, g_);
		boost::remove_vertex(v, g_);
	}
}

void UnitigGraph::removeStableSets()
{
	boost::graph_traits<UGraph>::vertex_iterator vi, vi_end, next;
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
void UnitigGraph::markCycles() //non-recusrive tarjan implementation for unitig graph
{
	auto cc = boost::get(boost::vertex_index2_t(), g_); // scc
	auto i = boost::get(boost::vertex_discover_time_t(), g_); // index
	auto onStack = boost::get(boost::vertex_finish_time_t(), g_); // onStack
    // lca, curr, scc
	std::stack<std::pair<std::pair<UVertex*, UVertex*>, unsigned int> > recursion_stack;
	std::stack<UVertex*> visit_stack;
	unsigned int index = 1;
	for (auto& v : boost::vertices(g_))
	{
        if (boost::get(cc,v) == 0) // scc hasnt been set
		{
			recursion_stack.push(std::make_pair(std::make_pair(nullptr,&v),0));
		}
		while (recursion_stack.size() > 0)
		{
			auto next_element = recursion_stack.top();
			recursion_stack.pop();

			UVertex *prev = next_element.first.first; // where we came from
			UVertex curr = *(next_element.first.second);
			
			unsigned int child = next_element.second;
			unsigned int children = boost::out_degree(curr,g_);//check reverse complement?
		    
			if (!child) // this vertex is visited the first time this run (child is nullptr)
			{
				boost::put(i,curr,index);
				boost::put(cc,curr,index);
				index++;
				boost::put(onStack,curr,true);
				visit_stack.push(&curr);
			}
			if (child < children) // still have to search at least one child
			{
                // reverse complements? TODO
                auto successors = boost::out_edges(curr,g_);
                unsigned int j = 0;
                UVertex next;
                for (auto&& ne : successors)
                {
                    j++;
                    next = boost::target(ne,g_);
                    if (j > child)
                        break;
                }
                ++child;
				recursion_stack.push(std::make_pair(std::make_pair(prev,&curr),child)); // visit the next child
				if (boost::get(i,next) == 0)
				{
					recursion_stack.push(std::make_pair(std::make_pair(&curr,&next),0)); // "recurse" on current child
				}
				else if (boost::get(onStack, next))
				{
                    boost::put(cc,curr,std::min(boost::get(cc,curr),boost::get(i,next)));
				}
			}
			else  // all children have been visited in the stack
			{
				if (prev) // is not the first searched, prev is not nullptr
				{
					boost::put(cc,*prev,std::min(boost::get(cc,*prev),boost::get(cc,curr)));
					if (boost::get(i,*prev) == boost::get(cc,*prev))
					{
						UVertex* scc = visit_stack.top();
						do
						{
							visit_stack.pop();
							boost::put(onStack,*scc,false);
							scc = visit_stack.top();
						} while (*scc != *prev);
					}
				}
				else
				{
					UVertex *top = visit_stack.top();
					boost::put(onStack,*top,false);
					visit_stack.pop();
				}
			}
		}
	}
}

// the graph might contain some unconnected vertices, clean up
void UnitigGraph::cleanGraph()
{
	std::cerr << "Removing stable sets..." << std::endl;
	removeStableSets();
	if (true) //TODO make this parametrized if results isnt always strictly improved
	{
		std::cerr << "Filtering terminal paths..." << std::endl;
		filterTerminals(); // this can be done in cycles (filter - contract - remove)
		std::cerr << "Contracting simple paths.." << std::endl;
		contractPaths(); 
		std::cerr << "Removing stable sets (again)..." << std::endl;
		removeStableSets();
	}
}

// returns the sources of the graph sorted by their connected components
std::vector<Connected_Component> UnitigGraph::getSources() const
{
	auto cc = boost::get(boost::vertex_index2_t(), g_);
	unsigned int cc_count = 0; // "real" number of ccs
	std::vector<std::vector<UVertex> > sources;
	std::unordered_map<unsigned int,unsigned int> cc_map;
	for (auto&& v : boost::vertices(g_))
	{
		if (boost::in_degree(v,g_) == 0)
		{
			unsigned int curr_cc = boost::get(cc,v);
			auto miter = cc_map.find(curr_cc);
			if (miter == cc_map.end())
			{
				std::vector<UVertex> cc_i({v});
				sources.push_back(cc_i);
				cc_map.insert(std::make_pair(curr_cc,cc_count++));
			}
			else
			{
				sources[miter->second].push_back(v);
			}
		}
	}
	return sources;
}

// add edges of source to heap
void UnitigGraph::add_sorted_edges(std::vector<UEdge>& q, const UVertex& source, bool addAll)
{
	const auto& cap = boost::get(boost::edge_capacity_t(),g_);
	
	auto edgeCompare = [&](UEdge e1, UEdge e2){
		return boost::get(cap, e1) < boost::get(cap, e2); // highest capacity first
	}; //lambda for comparing to edges based on their capacity
	
	auto out_edges = boost::out_edges(source, g_);
	for (const auto& oe : out_edges) // we dont do coverage checks here TODO?
	{
		q.push_back(oe);
		std::push_heap(q.begin(), q.end(), edgeCompare);
	}
}

// Tests whether two percentages "belong together" TODO this is quite arbitrary
bool UnitigGraph::test_hypothesis(float to_test, float h0)
{
	return (std::abs(to_test - h0) < FRAC_THRESH); // they differ by less than 6%
	// this probably is not bad for values significantly larger than 5%
}

void UnitigGraph::find_fattest_path(UVertex target, std::string& sequence, std::vector<float>& coverage_fraction, std::vector<UEdge>& visited_edges)
{
	const auto& cap = boost::get(boost::edge_capacity_t(),g_);
	const auto& ename = boost::get(boost::edge_name_t(), g_);
	const auto& visit = boost::get(boost::edge_index_t(), g_);
	
	int outdegree = boost::out_degree(target, g_);
	UEdge next_edge;

	while (outdegree) // not a sink
	{
		std::vector<UEdge> fattest_edges;
		std::make_heap(fattest_edges.begin(), fattest_edges.end());
		
		add_sorted_edges(fattest_edges, target, true);
		
		std::sort_heap(fattest_edges.begin(), fattest_edges.end());

		float total_coverage = 0.;
		
		for (const auto& e : fattest_edges)
		{
			total_coverage += boost::get(cap, e);
		}
		next_edge = fattest_edges.back();
		// if next_edge is visited, we found a cycle, break it by following the second fattest edge
		if (boost::get(visit, next_edge))
		{
			break;
			while (!fattest_edges.empty() and boost::get(visit, next_edge))
			{
				next_edge = fattest_edges.back();
				std::pop_heap(fattest_edges.begin(), fattest_edges.end());
				fattest_edges.pop_back();
			}
		}
		else
		{
			std::pop_heap(fattest_edges.begin(), fattest_edges.end());
			fattest_edges.pop_back();
		}
		target = boost::target(next_edge, g_);
		sequence += boost::get(ename, next_edge); // add to current sequence
		coverage_fraction.push_back(boost::get(cap,next_edge)/total_coverage); // fraction of outflow
		visited_edges.push_back(next_edge); // has already been updated here
		boost::put(visit, next_edge, true); // mark edge as visited
		outdegree = boost::out_degree(target, g_);
	}
	float min_fraction = *std::min_element(coverage_fraction.begin(), coverage_fraction.end());
	for (unsigned int i = 0; i < coverage_fraction.size(); i++)
	{
		if (!test_hypothesis(coverage_fraction[i],min_fraction)) // significantly higher, some flow remains
		{
			float new_cov = boost::get(cap, visited_edges[i]) * (1 - min_fraction);
			boost::put(cap, visited_edges[i], new_cov); // reduce by min_fraction percent
			boost::put(visit, visited_edges[i], false); // all edges with reamining coverage > 0 might be visited again
		}
		else // they belong to the same genome, use all flow
		{
			boost::put(cap, visited_edges[i], 0.);
			boost::put(visit, visited_edges[i], false); // "unvisit" edges
		}
	}
}

// calculates the flows and corresponding paths through the graph
void UnitigGraph::calculateFlow()
{
	std::cerr << "Calculating flow..." << std::endl;
	const auto& cap = boost::get(boost::edge_capacity_t(),g_);
	const auto& vname = boost::get(boost::vertex_name_t(), g_);
	const auto& ename = boost::get(boost::edge_name_t(), g_);

	// sources is a vector of all sources of a certain connected component
	auto sources = getSources();

	std::cerr << sources.size() << " connected components found." << std::endl;

	unsigned int j = 0; // contig number
	for (const Connected_Component& cc : sources)
	{
		std::vector<UEdge> q;
		std::make_heap(q.begin(), q.end());
		// sort all first edges by their capacity
		for (const auto& source : cc)
		{
			add_sorted_edges(q, source, false); //creates a heap
		}
		std::sort_heap(q.begin(), q.end());

		UEdge first_edge;
		while (!q.empty())
		{
			first_edge = q.back(); //

			float first_capacity = boost::get(cap, first_edge);

			std::vector<float> coverage_fraction; // fraction of the path
			std::vector<UEdge> visited_edges; // edges on the fattest path, TODO store pointers?

			auto source = boost::source(first_edge, g_);
			auto target = boost::target(first_edge, g_);
			std::string sequence = boost::get(vname, source) + boost::get(ename,first_edge);
			float total_coverage = 0.;
			for (const auto& e : boost::out_edges(source, g_))
			{
				total_coverage += boost::get(cap,e);
			}
			coverage_fraction.push_back(boost::get(cap,first_edge)/total_coverage);
			visited_edges.push_back(first_edge);
		
			find_fattest_path(target, sequence, coverage_fraction, visited_edges);

			if (sequence.length() > CONTIG_THRESH) // i.e. readsize
			{
				std::cout << ">Contig_" << j++ << " (" << *std::min_element(coverage_fraction.begin(), coverage_fraction.end()) << " of " << first_capacity << ")" << std::endl;
				std::cout << sequence << std::endl;
			}
			
			if (boost::get(cap, first_edge) == 0)
			{
				q.pop_back(); std::pop_heap(q.begin(), q.end()); // source flow has been depleted
			}
			std::sort_heap(q.begin(), q.end()); // capacity of first_edge has changed
		}
	}
}

void UnitigGraph::debug()
{
    //markCycles();

	// DEBUG
	auto&& numV = boost::num_vertices(g_);
	auto&& numE = boost::num_edges(g_);
	
	std::cerr << "Unitig graph has " << numV << " vertices and " << numE << " edges, starting cleaning" << std::endl;
	clock_t t = clock();

	//TODO move this out
	cleanGraph();
	
	numV = boost::num_vertices(g_);
	numE = boost::num_edges(g_);
	std::cerr << "Finished cleaning: Unitig graph has " << numV << " vertices and " << numE << " edges" << std::endl;
	std::cerr << "Cleaning took " << (clock() - t)/1000000. << " seconds." << std::endl;

    t = clock();
    std::cerr << "Finding Strongly Connected Components" << std::endl;
    markCycles();
    std::cerr << "Done finding sccs after " << (clock() -t)/100000. << " seconds" << std::endl;

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
			total_out += boost::get(cap,e.first);
		}

		for (auto&& ie : boost::in_edges(*vi,g_))
		{
			auto&& e = boost::edge(boost::source(ie,g_),*vi,g_);
			total_in += boost::get(cap,e.first);
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

    //boost::write_graphviz(std::cout, g_, boost::make_label_writer(boost::get(boost::vertex_index1_t(),g_)), boost::make_label_writer(boost::get(boost::edge_capacity_t(),g_)), boost::default_writer(), propmapIndex);
	boost::write_graphviz(std::cout, g_, boost::make_label_writer(boost::get(boost::vertex_name_t(),g_)), boost::make_label_writer(boost::get(boost::edge_name_t(),g_)), boost::default_writer(), propmapIndex);
	//boost::write_graphviz(std::cout, g_, boost::make_label_writer(boost::get(boost::vertex_index2_t(),g_)), boost::make_label_writer(boost::get(boost::edge_capacity_t(),g_)), boost::default_writer(), propmapIndex);
}
