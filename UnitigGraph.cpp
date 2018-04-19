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

// unitig graph for debugging purposes
UnitigGraph::UnitigGraph() : cc_(1), threshold_(0)
{
    std::vector<UVertex> vertices;
    for (unsigned int i = 1; i < 7; i++)
    {
        unsigned int index = 1;
        vertices.push_back(addVertex(&index, std::to_string(i)));
    }
    boost::add_edge(vertices[0],vertices[1],g_);
    boost::add_edge(vertices[0],vertices[2],g_);
    boost::add_edge(vertices[1],vertices[3],g_);
    boost::add_edge(vertices[1],vertices[4],g_);
    boost::add_edge(vertices[2],vertices[3],g_);
    boost::add_edge(vertices[2],vertices[4],g_);
    boost::add_edge(vertices[3],vertices[5],g_);
    boost::add_edge(vertices[4],vertices[5],g_);
    /*boost::add_edge(vertices[0],vertices[1],g_);
    boost::add_edge(vertices[0],vertices[1],g_);
    boost::add_edge(vertices[1],vertices[2],g_);
    boost::add_edge(vertices[1],vertices[3],g_);
    boost::add_edge(vertices[1],vertices[12],g_);
    boost::add_edge(vertices[2],vertices[3],g_);
    boost::add_edge(vertices[2],vertices[4],g_);
    boost::add_edge(vertices[3],vertices[0],g_);
    boost::add_edge(vertices[4],vertices[5],g_);
    //boost::add_edge(vertices[4],vertices[5],g_);
    boost::add_edge(vertices[4],vertices[8],g_);
    boost::add_edge(vertices[5],vertices[6],g_);
    boost::add_edge(vertices[6],vertices[7],g_);
    boost::add_edge(vertices[7],vertices[6],g_);
    boost::add_edge(vertices[7],vertices[8],g_);
    boost::add_edge(vertices[8],vertices[9],g_);
    boost::add_edge(vertices[9],vertices[10],g_);
    boost::add_edge(vertices[9],vertices[11],g_);
    boost::add_edge(vertices[10],vertices[14],g_);
    boost::add_edge(vertices[10],vertices[14],g_);
    boost::add_edge(vertices[11],vertices[12],g_);
    boost::add_edge(vertices[11],vertices[13],g_);
    boost::add_edge(vertices[12],vertices[4],g_);
    boost::add_edge(vertices[14],vertices[13],g_);
    //boost::add_edge(vertices[14],vertices[15],g_);*/
}
// constructor of the so-called UnitigGraph
// unifies all simple paths in the deBruijnGraph to a single source->sink path
// all remaining nodes have either indegree != outdegree or indegree == outdegree > 1
UnitigGraph::UnitigGraph(deBruijnGraph& dbg, float error_rate) : cc_(1), threshold_(0)
{
	std::cerr << "deBruijnGraph has " << dbg.getSize() << " vertices" << std::endl;
	std::cerr << "Building unitig graph from deBruijn graph..." << std::endl;
    
	clock_t t = clock();
	auto&& junc = dbg.getJunctions();
	unsigned int index = 1;
	auto&& out_unbalanced = junc.first;
	auto&& in_unbalanced = junc.second;
    threshold_ = calculateThresholds(dbg, error_rate);
	// starting from the sources, we build the unitig graph
	for (auto& v : out_unbalanced)
	{
		std::string curr = v.get_kmer();
		Vertex* source = dbg.getVertex(curr);
		// make a guess whether we are relevant already 
		if (!source->is_flagged())
		{
			connectUnbalanced(source, &index, curr, dbg, error_rate); 
			cc_++; // new connected component found
		}
	}
	for (auto& v : in_unbalanced)
	{
		std::string curr = v.get_kmer();
		Vertex* source = dbg.getVertex(curr);
		if (!source->is_flagged())
		{
			connectUnbalanced(source, &index, curr, dbg, error_rate); 
			cc_++; // new connected component found
		}
	}
	std::cerr << "Unitig graph succesfully build in " << (clock() - t)/1000000. << " seconds." << std::endl;
}

// TODO
float UnitigGraph::calculateThresholds(const deBruijnGraph& dbg, float error_rate)
{
    auto&& cov_distr = dbg.coverageDistribution();
    std::vector<std::pair<float, float> > sorted_coverage;
    float total_edges = 0.;
    float avg_coverage = 0.;
    float median_coverage = 0.;
    float middle_coverage = 0.;
    for (const auto& cov : cov_distr)
    {
        sorted_coverage.push_back(cov);
        total_edges += cov.second; //counts all edges
        avg_coverage += cov.first * cov.second; // adds up the coverage values
    }
    std::sort(sorted_coverage.begin(), sorted_coverage.end());
    float used_edges = 0.;
    unsigned int i = 0;
    for (const auto& cov : sorted_coverage)
    {
        used_edges += cov.second;
        if (used_edges > total_edges/2 and !median_coverage)
        {
            median_coverage = cov.first;
        }
        if (i++ > sorted_coverage.size()/2 and !middle_coverage)
        {
            middle_coverage = cov.first;
        }
    }
    avg_coverage /= total_edges;
    return middle_coverage * error_rate;
}

// adds a vertex to the unitig graph: adds it to the boost graph, as well as to the mapping from index to vertex
UVertex UnitigGraph::addVertex(unsigned int* index, std::string name)
{
	UVertex uv = boost::add_vertex(g_);
	(*index)++;
    // set vertex properties
    g_[uv].name = name;
	g_[uv].index = *index;
    g_[uv].scc = 1;
    g_[uv].tarjan_index = 0; // needs to be 0 to find out whether it has been set
    g_[uv].onStack = false;
    g_[uv].cc = cc_; // set cc to current CC
	
    auto&& ins = std::make_pair(*index,uv);
	graph_.insert(ins);
	return uv;
}

// function for connecting a given source/sink vertex to all its unbalanced successors/predecessors
void UnitigGraph::connectUnbalanced(Vertex* source, unsigned int* index, std::string curr, deBruijnGraph& dbg, float error_rate)
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
            if (junction->get_out_coverage(c) > error_rate*junction->get_total_out_coverage() and junction->get_out_coverage(c) > threshold_)
			{
				to_search = true; break;
			}
		}
		if (!to_search)
		{
			for (const char& c : pred)
			{
				if (junction->get_in_coverage(c) > error_rate * junction->get_total_in_coverage() and junction->get_in_coverage(c) > threshold_)
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
	float epsilon = float(total_out)/total_in; // the gain/loss of this vertex
	Sequence src = *dbg.getSequence(curr);
	// this is for checking whether the reverse complement is in the debruijn graph
	bool reverse = (src != curr); // true, if reverse complement, i.e. src and curr are not the same
	if (total_out) // out_degree >= 1
    {
        for (const auto& n : succ)
        {
            std::string sequence("");
            std::string next;
            float coverage = currV->get_out_coverage(n);
            float pcov = coverage/float(total_out); // percentage of coverage
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
    }
	// finding the predecessing unbalanced vertices
    if (total_in)
    {
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
    float first = coverage; // first out-coverage
    float last = coverage;
	unsigned int length = 1;
	char lastchar = g_[trg].name.back(); //char with which we are pointing to trg
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
        last = cov;
		length++; // used for the coverage caluclations of the path later on
		nextV = dbg.getVertex(prev);
		pred = nextV->get_predecessors();
		succ = nextV->get_successors();
		sequence += lastchar;
	}
	avg /= float(length); // average coverage over the path
	if (!nextV->is_visited() and avg >= threshold_) // if the next vertex has been visited it already is part of the unitiggraph, otherwise add it
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
    else if (avg < threshold_)
    {
        return std::make_pair(nextV,"");
    }
	// if the vertex has not been added, dont try to add an edge to it
	if (nextV->index == 0)
	{
		return std::make_pair(nextV,"");
	}
	UVertex src = graph_[nextV->index];
	auto e = boost::edge(src,trg,g_);
	// if edge has been added or the immediate neighbour is an unbalanced vertex, do not add edge
	if ((!e.second or (e.second and g_[e.first].name != sequence)))
	{
		e = boost::add_edge(src,trg,g_);
		std::reverse(sequence.begin(), sequence.end()); // we add the path from the found node to trg, the sequence was added in reverse order
		std::string old_name = g_[e.first].name;
		if (e.second) // TODO v-S->w-T->v is treated like v<-S-w<-T-v (should be ST self-loop and TS self-loop)
		{
            g_[e.first].name = old_name + sequence;
		}
		else
		{
            g_[e.first].name = sequence;    
		}
        g_[e.first].capacity = avg;
        g_[e.first].cap_info.avg = avg;
        g_[e.first].cap_info.max = max;
        g_[e.first].cap_info.min = min;
        g_[e.first].cap_info.first = first;
        g_[e.first].cap_info.last= last;
        g_[e.first].cap_info.length = g_[e.first].name.length();
        g_[e.first].visited = false;
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
    float first = coverage;
    float last = coverage;
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
        last = cov;
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
	if (!nextV->is_visited() and avg >= threshold_)
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
        // TODO this should not happen
		return std::make_pair(nextV,""); // path has been found, do not add anything
	}
    else if (avg < threshold_)
    {
        return std::make_pair(nextV,"");
    }
	// this vertex has been found and considered irrelevant because of insufficient coverage
	if (nextV->index == 0)
	{
		return std::make_pair(nextV,"");
	}
	
	UVertex trg = graph_[nextV->index];
	auto e = boost::edge(src,trg,g_);
	if ((!e.second or (e.second and g_[e.first].name != sequence)))
	{
		e = boost::add_edge(src,trg,g_);
		auto old_name = g_[e.first].name;
		if (e.second) // TODO v-S->w-T->v is treated like v<-S-w<-T-v (should be ST self-loop and TS self-loop)
		{
            g_[e.first].name = sequence + old_name;
		}
		else
		{
            g_[e.first].name = sequence;
		}
        g_[e.first].capacity = avg;
        g_[e.first].cap_info.avg = avg;
        g_[e.first].cap_info.max = max;
        g_[e.first].cap_info.min = min;
        g_[e.first].cap_info.first = first;
        g_[e.first].cap_info.last = last;
        g_[e.first].cap_info.length = g_[e.first].name.length();
        g_[e.first].visited = false;
	}
	return std::make_pair(nextV,next);
}

// contracts all simple paths in graph to a single source-sink connection
void UnitigGraph::contractPaths()
{
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
            auto&& f = boost::edge(*vi, new_target,g_); // coverage etc of second edge to be contracted
			std::string seq = g_[e.first].name;
			unsigned int w = seq.length();

			Capacity cap_info_e = g_[e.first].cap_info;
			Capacity cap_info_f = g_[f.first].cap_info;
            float max = std::max(cap_info_e.max, cap_info_f.max);
            float min = std::min(cap_info_e.min, cap_info_f.min);
            float first = cap_info_e.first;
            float last = cap_info_f.last;
            
            float capacity = g_[e.first].capacity * w;
			e = boost::edge(*vi,new_target,g_);
			seq += g_[e.first].name; // append the sequence
			capacity += g_[e.first].capacity * (seq.length() - w);
			capacity /= seq.length(); // currently using the average coverage on the contracted path

			auto&& new_e = boost::add_edge(new_source,new_target,g_);
            g_[new_e.first].name = seq;
            g_[new_e.first].capacity = capacity;
            g_[new_e.first].cap_info.avg = capacity;
            g_[new_e.first].cap_info.max = max;
            g_[new_e.first].cap_info.min = min;
            g_[new_e.first].cap_info.first = first;
            g_[new_e.first].cap_info.last = last;
            g_[new_e.first].cap_info.length = g_[new_e.first].name.length();
            g_[new_e.first].visited = false;
			boost::clear_vertex(*vi,g_);
			boost::remove_vertex(*vi,g_);
		}
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
    // lca, curr, scc
	std::stack<std::pair<std::pair<UVertex, UVertex>, unsigned int> > recursion_stack;
	std::stack<UVertex> visit_stack;
    UVertex dummy = boost::add_vertex(g_); // dummy vertex for starting from
    g_[dummy].scc = 0; // cc for all real vertices is at least 1
 	unsigned int index = 1;
	for (auto& v : boost::vertices(g_))
	{
        if (g_[v].tarjan_index == 0) // scc hasnt been set
		{
			recursion_stack.push(std::make_pair(std::make_pair(dummy,v),0));
		}
		while (recursion_stack.size() > 0)
		{
			auto next_element = recursion_stack.top();
			recursion_stack.pop();

			UVertex prev = next_element.first.first; // where we came from
			UVertex curr = next_element.first.second;

			unsigned int child = next_element.second;
			unsigned int children = boost::out_degree(curr,g_); //check reverse complement?
			
            if (!child) // this vertex is visited the first time this run (child is 0)
			{
                g_[curr].tarjan_index = index;
                g_[curr].scc = index;
                g_[curr].onStack = true;
				index++;
				visit_stack.push(curr);
			}
			if (child < children) // still have to search at least one child
			{
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
                
				recursion_stack.push(std::make_pair(std::make_pair(prev,curr),child)); // visit the next child
				if (g_[next].tarjan_index == 0)
				{
					recursion_stack.push(std::make_pair(std::make_pair(curr,next),0)); // "recurse" on current child
				}
				else if (g_[next].onStack)
				{
                    //g_[curr].scc = std::min(g_[curr].scc,g_[next].tarjan_index); // doing this seems to find SCCs within SCCs - might also be useful
                    g_[curr].scc = std::min(g_[curr].scc,g_[next].scc);
				}
			}
			else  // all children have been visited in the stack
			{
				if (g_[curr].scc) // is not the first searched
				{
					g_[prev].scc = std::min(g_[prev].scc,g_[curr].scc);
					if (g_[curr].scc == g_[curr].tarjan_index)
					{
						UVertex scc = visit_stack.top();
                        g_[scc].onStack = false;
                        visit_stack.pop();
                        while (g_[curr].name != g_[scc].name) // there should not be two vertices with the same name present
						{
							scc = visit_stack.top();
                            g_[scc].onStack = false;
                            visit_stack.pop();
						} 
					}
				}
				else
				{
					UVertex top = visit_stack.top();
                    g_[top].onStack = false;
					visit_stack.pop();
				}
			}
		}
	}
    boost::remove_vertex(dummy, g_); // remove the dummy again
}

// the graph might contain some unconnected vertices, clean up
void UnitigGraph::cleanGraph()
{
	std::cerr << "Removing stable sets..." << std::endl;
	removeStableSets();
    std::cerr << "Contracting simple paths.." << std::endl;
    contractPaths(); 
}

// returns the sources of the graph sorted by their connected components
std::vector<Connected_Component> UnitigGraph::getSources() const
{
	unsigned int cc_count = 0; // "real" number of ccs
	std::vector<std::vector<UVertex> > sources;
	std::unordered_map<unsigned int,unsigned int> cc_map;
	for (auto&& v : boost::vertices(g_))
	{
		if (boost::in_degree(v,g_) == 0)
		{
			unsigned int curr_cc = g_[v].cc;
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

// Tests whether two percentages "belong together" TODO this is quite arbitrary
bool UnitigGraph::test_hypothesis(float to_test, float h0)
{
	return (std::abs(to_test - h0) < 0.06); // they differ by less than 6%
    // TODO
}

void UnitigGraph::find_fattest_path(UVertex target, std::string& sequence, std::vector<std::pair<float, float> >& coverage_fraction, std::vector<UEdge>& visited_edges)
{
    auto edgeCompare = [&](UEdge e1, UEdge e2){
        return g_[e1].capacity < g_[e2].capacity; // highest capacity last
    }; //lambda for comparing to edges based on their capacity
	
    int outdegree = boost::out_degree(target, g_);
    std::vector<UEdge> cycles;
    UEdge fattest_edge;
    float fattest_coverage = 0.;
    bool cycle = false; // are we currently on a cycle
	float min_fraction = 1.01; //slightly above 1 so it is true for 1.0
    float max_capacity = 0.;
	UEdge next_edge;

	while (outdegree) // not a sink
	{
        UVertex curr = target;
		std::vector<UEdge> fattest_edges;
        auto out_edges = boost::out_edges(target, g_);
        for (const auto& oe : out_edges)
        {
            fattest_edges.push_back(oe);
        }
		std::sort(fattest_edges.begin(), fattest_edges.end(), edgeCompare);

		float total_coverage = 0.;
		
		for (const auto& e : fattest_edges)
		{
			total_coverage += g_[e].capacity;
		}
		next_edge = fattest_edges.back();
		target = boost::target(next_edge, g_);
		// if next_edge is visited, we found a cycle, "break" it
		if (g_[target].scc == g_[curr].scc) // a cycle starts here
            cycle = true;
        else //fattest edge does not pass through cycle
        {
            cycle = false;
            cycles.clear();
        }
        if (g_[next_edge].visited) //fattest path goes along cycle - now perform magic
		{
		}
		else
		{
			fattest_edges.pop_back();
		}
        if (cycle)
        {
			while (!fattest_edges.empty())
			{
				auto&& alt_edge = fattest_edges.back();
				auto&& alt_vertex = boost::target(alt_edge, g_);
				fattest_edges.pop_back();
                if (g_[alt_vertex].scc != g_[target].scc)
                {
                    cycles.push_back(alt_edge); // add all non-cycle edges
                    if (g_[alt_edge].capacity > fattest_coverage)
                    {
                        fattest_coverage = g_[alt_edge].capacity;
                        fattest_edge = alt_edge; // this will be the edge with which we are leaving the cycle
                    }
                    break;
                }
			}
        }
        float capacity = g_[next_edge].capacity;
		float fraction = capacity/total_coverage;
        if (test_hypothesis(fraction, 0.5)) // fraction is close to 0.5 -> split contigs?
        {
            // TODO
        }
		sequence += g_[next_edge].name; // add to current sequence
        coverage_fraction.push_back(std::make_pair(capacity, total_coverage)); // fraction of outflow
        if (fraction < min_fraction)
            min_fraction = fraction;
        if (capacity > max_capacity)
            max_capacity = capacity;
		visited_edges.push_back(next_edge); // has already been updated here
        g_[next_edge].visited = true; // mark edge as visited
		outdegree = boost::out_degree(target, g_);
	}
	for (unsigned int i = 0; i < coverage_fraction.size(); i++)
	{
		if (!test_hypothesis(coverage_fraction[i].first/coverage_fraction[i].second,min_fraction)) // significantly higher, some flow remains
		{
			float new_cov = g_[visited_edges[i]].capacity * (1 - min_fraction);
            g_[visited_edges[i]].capacity = new_cov; // reduce by min_fraction percent
            g_[visited_edges[i]].visited = false; // all edges with reamining coverage > 0 might be visited again
		}
		else // they belong to the same genome, use all flow
		{
            g_[visited_edges[i]].capacity = 0.;
            g_[visited_edges[i]].visited = false; // "unvisit" edges
		}
	}
    if (sequence.length() > 150) // i.e. readsize TODO
    {
        std::cout << ">Contig_(" << min_fraction << " of " << max_capacity << ")" << std::endl;
        std::cout << sequence << std::endl;
    }
}
	

// calculates the flows and corresponding paths through the graph
void UnitigGraph::calculateFlow()
{
    auto edgeCompare = [&](UEdge e1, UEdge e2){
        return g_[e1].capacity < g_[e2].capacity; // highest capacity last
    }; //lambda for comparing to edges based on their capacity
	
    std::cerr << "Calculating flow..." << std::endl;
	// sources is a vector of all sources of a certain connected component
	auto sources = getSources();

	for (const Connected_Component& cc : sources)
	{
		std::vector<UEdge> q;
		// sort all first edges by their capacity
		for (const auto& source : cc)
		{
            auto out_edges = boost::out_edges(source, g_);
            for (const auto& oe : out_edges)
            {
                q.push_back(oe);
		    }
        }
		std::sort(q.begin(), q.end(), edgeCompare);

		UEdge first_edge;
		while (!q.empty())
		{
			first_edge = q.back(); //

			float first_capacity = g_[first_edge].capacity;
            if (first_capacity < threshold_) // the highest abundant source is not relevant
            {
                q.pop_back();
                continue;
                //break;
            }

			std::vector<std::pair<float, float> > coverage_fraction; // fraction of the path
			std::vector<UEdge> visited_edges; // edges on the fattest path, TODO store pointers?

			auto source = boost::source(first_edge, g_);
			auto target = boost::target(first_edge, g_);
			std::string sequence = g_[source].name + g_[first_edge].name;
			float total_coverage = 0.;
			for (const auto& e : boost::out_edges(source, g_))
			{
				total_coverage += g_[e].capacity;
			}
			coverage_fraction.push_back(std::make_pair(g_[first_edge].capacity,total_coverage));
			visited_edges.push_back(first_edge);
		
			find_fattest_path(target, sequence, coverage_fraction, visited_edges);
			
			if (g_[first_edge].capacity == 0)
			{
				q.pop_back(); // source flow has been depleted
			}
            /* for statistical evaluation of the resulting pathes*/
			std::sort(q.begin(), q.end(), edgeCompare); // capacity of first_edge has changed
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
    std::cerr << "Done finding sccs after " << (clock() - t)/100000. << " seconds" << std::endl;

	typedef std::map<UVertex, int> IndexMap;
	IndexMap mapIndex;
	boost::associative_property_map<IndexMap> propmapIndex(mapIndex);
	uvertex_iter vi, vi_end;
	int i = 1;
	int simpletons = 0;
	for (boost::tie(vi,vi_end) = boost::vertices(g_); vi != vi_end; ++vi)
	{
		boost::put(propmapIndex,*vi,i++);
	}

    //boost::write_graphviz(std::cout, g_, boost::make_label_writer(boost::get(&VertexProperties::index,g_)), boost::make_label_writer(boost::get(&EdgeProperties::capacity,g_)), boost::default_writer(), propmapIndex);
	//boost::write_graphviz(std::cout, g_, boost::make_label_writer(boost::get(&VertexProperties::name, g_)), boost::make_label_writer(boost::get(&EdgeProperties::name,g_)), boost::default_writer(), propmapIndex);
	//boost::write_graphviz(std::cout, g_, boost::make_label_writer(boost::get(&VertexProperties::cc,g_)), boost::make_label_writer(boost::get(&EdgeProperties::capacity,g_)), boost::default_writer(), propmapIndex);
}
