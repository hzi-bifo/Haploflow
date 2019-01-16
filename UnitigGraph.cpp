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
    threshold_ = calculate_thresholds(dbg, error_rate);
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
	std::cerr << "Unitig graph successfully build in " << (clock() - t)/1000000. << " seconds." << std::endl;
    t = clock();
    std::cerr << "Finding Strongly Connected Components" << std::endl;
    markCycles();
    std::cerr << "Done finding sccs after " << (clock() - t)/100000. << " seconds" << std::endl;
}

float UnitigGraph::calculate_thresholds(const deBruijnGraph& dbg, float error_rate)
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
    return middle_coverage * error_rate; //TODO
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
    g_[uv].visiting_time = 0;
	
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
                float curr_start = currV->get_read_starts() * currV->get_out_coverage(n)/currV->get_total_out_coverage();
                next = curr.substr(1) + n;
                Vertex* nextV = dbg.getVertex(next);
                Sequence s = *dbg.getSequence(next);
                sequence += n;
                if (!nextV->is_flagged())
                    following.push_back(buildEdge(uv, nextV, next, sequence, index, coverage, pcov, dbg, curr_start));
            }
            // if we are a reverse complement, we actually want to add the path in reverse order
            else
            {
                float curr_end = currV->get_read_ends() * currV->get_in_coverage(deBruijnGraph::complement(n))/currV->get_total_in_coverage();
                next = deBruijnGraph::complement(n) + curr.substr(0,curr.length() - 1);
                Vertex* nextV = dbg.getVertex(next);
                Sequence s = *dbg.getSequence(next);
                sequence += curr.back(); // the predecessor points to the current vertex with the last char of curr (by definition)
                if (!nextV->is_flagged())
                    following.push_back(buildEdgeReverse(uv, nextV, next, sequence, index, coverage, pcov, dbg, curr_end));
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
                float curr_end = currV->get_read_ends() * currV->get_in_coverage(n)/currV->get_total_in_coverage();
                prev = n + curr.substr(0,curr.length() - 1);
                Vertex* nextV = dbg.getVertex(prev);
                Sequence s = *dbg.getSequence(prev);
                sequence += curr.back(); // the predecessor points to the current vertex with the last char of curr
                if (!nextV->is_flagged())
                    following.push_back(buildEdgeReverse(uv, nextV, prev, sequence, index, coverage, pcov, dbg, curr_end));
            }
            else
            {
                float curr_start = currV->get_read_starts() * currV->get_out_coverage(deBruijnGraph::complement(n))/currV->get_total_out_coverage();
                prev = curr.substr(1) + deBruijnGraph::complement(n);
                Vertex* nextV = dbg.getVertex(prev);
                Sequence s = *dbg.getSequence(prev);
                sequence += deBruijnGraph::complement(n);
                if (!nextV->is_flagged())
                    following.push_back(buildEdge(uv, nextV, prev, sequence, index, coverage, pcov, dbg, curr_start));
            }
        }
    }
	currV->flag(); // this vertex is done
	return following;
}

// go back through the graph until the next unbalanced node is found and add an ("reversed") edge
std::pair<Vertex*,std::string> UnitigGraph::buildEdgeReverse(UVertex trg, Vertex* nextV, std::string prev, std::string& sequence, unsigned int* index, float coverage, float pcov, deBruijnGraph& dbg, float curr_end)
{
	char lastchar = g_[trg].name.back(); //char with which we are pointing to trg
    float total_out = nextV->get_total_out_coverage();
    if (total_out == 0)
    {
        unsigned int comp = 1;
        total_out = std::max(nextV->get_out_coverage(lastchar),comp);
    }
    float starts_with = nextV->get_read_starts() * nextV->get_out_coverage(lastchar)/total_out; // number of reads starts within this edge, normalized by total flow through this edge
    float ends_with = curr_end; // number of reads ends within this edge
	auto&& succ = nextV->get_successors();
	auto&& pred = nextV->get_predecessors();
	// DEBUG
	float min = coverage; // TODO
	float max = coverage;
	float avg = coverage;
    float first = coverage; // first out-coverage
    float last = coverage;
	unsigned int length = 1;
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
        first = cov;
		length++; // used for the coverage caluclations of the path later on
		nextV = dbg.getVertex(prev);
		pred = nextV->get_predecessors();
		succ = nextV->get_successors();
		sequence += lastchar;
        starts_with += nextV->get_read_starts();
        ends_with += nextV->get_read_ends();
        if (std::abs(last - cov) > threshold_)
            break;
	}
	avg /= float(length); // average coverage over the path
	if (!nextV->is_visited() and avg >= threshold_) // if the next vertex has been visited it already is part of the unitiggraph, otherwise add it
	{
		nextV->visit();
		addVertex(index, prev); // the vertex is new and found to be relevant
		nextV->index = *index;
	}
	else if (!nextV->is_visited() or avg < threshold_ or nextV->index == 0)
	{
		return std::make_pair(nextV,""); // path has too low coverage
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
        g_[e.first].last_visit = 0;
        g_[e.first].capacity = avg;
        g_[e.first].cap_info.avg = avg;
        g_[e.first].cap_info.max = max;
        g_[e.first].cap_info.min = min;
        g_[e.first].cap_info.first = first;
        g_[e.first].cap_info.last = last;
        g_[e.first].cap_info.length = g_[e.first].name.length();
        g_[e.first].cap_info.starting = starts_with/(g_[e.first].cap_info.length * avg);
        g_[e.first].cap_info.ending = ends_with/(g_[e.first].cap_info.length * avg);
        g_[e.first].visited = false;
	}
	return std::make_pair(nextV,prev);
}

// same function like the reverse one, but going forward and finding successors
std::pair<Vertex*,std::string> UnitigGraph::buildEdge(UVertex src, Vertex* nextV, std::string next, std::string& sequence, unsigned int* index, float coverage, float pcov, deBruijnGraph& dbg, float curr_start)
{
	// with a little effort this can be moved inside the while loop for efficiency reasons
    float first_char = g_[src].name.front();
    float starts_with = curr_start; // number of reads starts within this edge, normalized by total flow through this edge
    float total_in = nextV->get_total_in_coverage();
    if (total_in == 0)
    {
        unsigned int comp = 1;
        total_in = std::max(nextV->get_in_coverage(first_char), comp);
    }
    float ends_with = nextV->get_read_ends() * nextV->get_in_coverage(first_char)/total_in; // number of reads ends within this edge
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
        starts_with += nextV->get_read_starts();
        ends_with += nextV->get_read_ends();
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
	else if (!nextV->is_visited() or avg < threshold_ or nextV->index == 0)
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
        g_[e.first].last_visit = 0;
        g_[e.first].capacity = avg;
        g_[e.first].cap_info.avg = avg;
        g_[e.first].cap_info.max = max;
        g_[e.first].cap_info.min = min;
        g_[e.first].cap_info.first = first;
        g_[e.first].cap_info.last = last;
        g_[e.first].cap_info.length = g_[e.first].name.length();
        g_[e.first].cap_info.starting = starts_with/(g_[e.first].cap_info.length * avg);
        g_[e.first].cap_info.ending = ends_with/(g_[e.first].cap_info.length * avg);
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
		
		// if in and outdegree is 1, we are on a simple path and can contract again
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

            // TODO temporarily disabled to study effect of not contracting
            //if (std::abs(cap_info_f.first - cap_info_e.last) > threshold_ or std::abs(cap_info_f.first - cap_info_e.last) > threshold_) // do not contract paths which have high divergence in capacity
            //    continue;
            
            float capacity = g_[e.first].capacity * w;
			e = boost::edge(*vi,new_target,g_);
			seq += g_[e.first].name; // append the sequence
			capacity += g_[e.first].capacity * (seq.length() - w);
			capacity /= seq.length(); // currently using the average coverage on the contracted path

			auto&& new_e = boost::add_edge(new_source,new_target,g_);
            g_[new_e.first].last_visit = 0;
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
void UnitigGraph::markCycles() //non-recursive tarjan implementation for unitig graph
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

void UnitigGraph::removeEmpty()
{
    std::vector<UEdge> toDelete;
    for (auto&& e : boost::edges(g_))
    {
        if (g_[e].capacity == 0)
            toDelete.push_back(e);
    }
    for (auto&& e : toDelete)
    {
        boost::remove_edge(e, g_);
    }
}

// the graph might contain some unconnected vertices, clean up
void UnitigGraph::cleanGraph()
{
    std::cerr << "Removing edges with no capacity..." << std::endl;
	removeEmpty();
    std::cerr << "Removing stable sets..." << std::endl;
	removeStableSets();
    //std::cerr << "Contracting simple paths..." << std::endl;
    //contractPaths(); 
}

// returns the seed vertex = vertex with the highest "coverage"
UEdge UnitigGraph::getSeed() const
{
    UEdge seed;
    UEdge source1;
    UEdge source2;
    UEdge source3;
    float max = 0;
    float max_diff = 0;
    float max_rel = 0;
    float max_start = 0;
	for (auto&& e : boost::edges(g_))
	{
        auto cap = g_[e].cap_info;
        float strt = cap.starting;
        float abs = std::abs(cap.first - cap.last);
        float rel = std::max(cap.first, cap.last) / std::min(cap.first, cap.last);
        if (abs > max_diff)
        {
            max_diff = abs;
            source1 = e;
        }
        if (rel > max_rel)
        {
            max_rel = rel;
            source2 = e;
        }
        if (strt > max_start)
        {
            max_start = strt;
            source3 = e;
        }
        if (g_[e].capacity > max)
        {
            max = g_[e].capacity;
            seed = e;
        }
    }
    auto src1 = boost::source(source1, g_);
    auto src2 = boost::source(source2, g_);
    auto src3 = boost::source(source3, g_);
    //std::cerr << g_[src1].index << " " << max_diff << std::endl;
    //std::cerr << g_[src2].index << " " << max_rel << std::endl;
    //std::cerr << g_[src3].index << " " << max_start << std::endl;
	return seed;
}

// Tests whether two percentages "belong together"
bool UnitigGraph::test_hypothesis(float to_test_num, float to_test_denom, float h0)
{
    float diff = std::abs(to_test_num - h0 * to_test_denom); // the absolute difference between num and expected num
    //float perc = to_test_num/to_test_denom; //currently unused
    return (diff < threshold_); // they differ by less than threshold 
}


float UnitigGraph::in_capacity(UVertex source)
{
    float capacity = 0;
    for (auto ie : boost::in_edges(source, g_))
    {
        capacity += g_[ie].capacity;
    }
    return capacity;
}

float UnitigGraph::out_capacity(UVertex target)
{
    float capacity = 0;
    for (auto oe : boost::out_edges(target, g_))
    {
        capacity += g_[oe].capacity;
    }
    return capacity;
}

std::pair<float, std::vector<float> > UnitigGraph::calculate_flow(std::vector<UEdge>& path)
{
    float min_flow = -1;
    std::vector<float> flows;

    unsigned int tot_length = 0;
    for (auto& e : path)
    {
        tot_length += g_[e].name.size(); // to see if calculated value is more than 150 bases away from begin
    }
    unsigned int length = 0;
    for (auto& e : path)
    {
        if ((length >= 150 and length <= tot_length - 150))
        {
            if (min_flow == -1 or g_[e].capacity < min_flow)
            {
                min_flow = g_[e].capacity;
            }
        }
        length += g_[e].name.size();
    }
    for (auto& e : path)
    {
        if (min_flow == -1) // TODO removes all flow (?)
        {
            flows.push_back(g_[e].capacity);
            continue;
        }
        flows.push_back(min_flow);
    } 
    return std::make_pair(1.0, flows);
}

// run dijsktra with fatness as optimality criterion
void UnitigGraph::dijkstra(UEdge seed, bool forward)
{
    auto edge_compare = [&](UEdge e1, UEdge e2){ //sort by biggest fatness
        if (g_[e1].fatness == g_[e2].fatness)
            return g_[e1].fatness2 < g_[e2].fatness2;
        else
            return g_[e1].fatness < g_[e2].fatness;
    };
    auto src = boost::source(seed, g_);
    auto trg = boost::target(seed, g_);
    std::vector<UEdge> q;
    for (auto e : boost::edges(g_)) //initialise distances and fatness
    {
        if (forward)
        {
            g_[e].fatness = 0;
            g_[e].distance = std::numeric_limits<unsigned int>::max();
        }
        else
        {
            g_[e].fatness2 = 0;
            g_[e].distance2 = std::numeric_limits<unsigned int>::max();
        }
        if (e == seed and forward)
        {
            g_[e].distance = 0;
            g_[e].fatness = g_[e].capacity;
        }
        else if (e == seed and !forward)
        {
            g_[e].distance2 = 0;
            g_[e].fatness2 = g_[e].capacity;
        }
    }
    q.push_back(seed);
    std::sort(q.begin(), q.end(), edge_compare);
    while (!q.empty()) // classic dijsktra routine
    {
        auto curr = q.back();
        q.pop_back();
        g_[curr].visited = true;
        auto source = boost::source(curr, g_);
        auto target = boost::target(curr, g_);
        if (forward)
        {
            for (auto oe : boost::out_edges(target, g_))
            {
                float fat = g_[oe].fatness;
                if (fat < std::min(g_[curr].fatness, g_[oe].capacity))
                {
                    g_[oe].fatness = std::min(g_[curr].fatness, g_[oe].capacity);
                    g_[oe].prev = curr;
                    g_[oe].distance = g_[curr].distance + g_[oe].name.size();
                }
                if (!g_[oe].visited)
                    q.push_back(oe);
            }
        }
        else
        {
            for (auto ie : boost::in_edges(source, g_))
            {
                float fat = g_[ie].fatness2;
                if (fat < std::min(g_[curr].fatness2, g_[ie].capacity))
                {
                    g_[ie].fatness2 = std::min(g_[curr].fatness2, g_[ie].capacity);
                    g_[ie].next = curr;
                    g_[ie].distance2 = g_[curr].distance2 + g_[ie].name.size();
                }
                if (!g_[ie].visited)
                    q.push_back(ie);
            }
        }
        std::sort(q.begin(), q.end(), edge_compare);
    }
    for (auto e : boost::edges(g_))
    {
        g_[e].visited = false;
    }
}

// Calculates the fattest path through the graph
std::vector<UEdge> UnitigGraph::find_fattest_path(UEdge seed)
{
    auto src = boost::source(seed, g_);
    auto trg = boost::target(seed, g_);
    std::cerr << g_[src].index << " --> " << g_[trg].index << " is seed" << std::endl;
    unvisit(); // make sure we start from unvisited graph
    dijkstra(seed, true);
    dijkstra(seed, false);
    float max_dist = 0;
    UEdge last;
    bool forward = true;
    for (auto e : boost::edges(g_))
    {
        auto source = boost::source(e, g_);
        auto target = boost::target(e, g_);
        if (g_[e].distance > max_dist and g_[e].distance < std::numeric_limits<unsigned int>::max())
        {
            max_dist = g_[e].distance; // so we have both maximised distance and fatness
            // TODO: magnitude of flow vs distance?
            last = e;
            forward = true;
        }
        if (g_[e].distance2 > max_dist and g_[e].distance2 < std::numeric_limits<unsigned int>::max())
        {
            max_dist = g_[e].distance2;
            last = e;
            forward = false;
        }
    }
    if (max_dist == 0) // path is only one edge, no longest path
    {
        return std::vector<UEdge>{seed};
    }
    auto curr = last;
    std::deque<UEdge> path = {curr};
    if (forward)
    {
        auto source = boost::source(curr, g_);
        while(g_[curr].distance > 0 and g_[curr].distance < std::numeric_limits<unsigned int>::max() and !g_[curr].visited) // this means the distance has been set, i.e. the vertex has been reached
        {
            g_[curr].visited = true;
            curr = g_[curr].prev;
            source = boost::source(curr, g_);
            path.push_front(curr);
        }
    }
    else
    {
        auto target = boost::target(curr, g_);
        while(g_[curr].distance2 > 0 and g_[curr].distance2 < std::numeric_limits<unsigned int>::max() and !g_[curr].visited) // this means the distance has been set, i.e. the vertex has been reached
        {
            g_[curr].visited = true;
            curr = g_[curr].next;
            target = boost::target(curr, g_);
            path.push_back(curr);
        }
    }
    src = boost::source(path.front(), g_);
    std::cerr << g_[src].index;
    for (auto& e : path)
    {
        auto target = boost::target(e, g_);
        std::cerr << " --> " << g_[target].index << " (" << std::max(g_[e].fatness, g_[e].fatness2) << ")";
    }
    std::cerr << std::endl;
    return std::vector<UEdge>(path.begin(), path.end());
}

void UnitigGraph::unvisit()
{
    for (auto& v : boost::vertices(g_))
    {
        g_[v].visiting_time = 0;
    }
    for (auto e : boost::edges(g_))
    {
        g_[e].visited = false;
        //g_[e].visits.visits = std::vector<unsigned int>{};
        //g_[e].last_visit = 0;
    }
}

// Calculates how much "gain" in flow a single vertex has in percent
float UnitigGraph::calculate_gain(UVertex& v)
{
    if (in_capacity(v) == 0 or out_capacity(v) == 0)
    {
        return 1; // no gain in source (definition)
    }
    auto inedges = boost::in_edges(v, g_);
    auto outedges = boost::out_edges(v, g_);
    float inflow = 0;
    float outflow = 0;
    for (auto&& ie : inedges)
    {
        inflow += g_[ie].cap_info.last;
    }
    for (auto&& oe : outedges)
    {
        outflow += g_[oe].cap_info.first;
    }
    return outflow/inflow;
}

// Given all the chosen edges and their coverage fraction, builds the contigs and reduces flow accordingly
std::pair<std::string, std::pair<float, float> > UnitigGraph::calculate_contigs(std::vector<UEdge>& path)
{
    auto flow = calculate_flow(path);
    float flow_frac = flow.first;
    std::vector<float> flows = flow.second;
    float min_flow = flows.front();

    unsigned int i = 0;
    UEdge start = path.front();
    UVertex source = boost::source(start, g_);
    std::string contig = g_[source].name; // kmer of first vertex
    for (auto& e : path)
    {   
        contig += g_[e].name;
        float current_flow = flows[i];
        if (test_hypothesis(g_[e].capacity, current_flow, 1.0) or g_[e].capacity < current_flow)
        {
            g_[e].capacity = 0;
            g_[e].cap_info.first = 0;
            g_[e].cap_info.last= 0;
            g_[e].cap_info.min = 0;
            g_[e].cap_info.max = 0;
            g_[e].cap_info.avg = 0;
        }
        else
        {
            g_[e].capacity -= current_flow;
            g_[e].cap_info.first -= current_flow;
            g_[e].cap_info.last -= current_flow;
            g_[e].cap_info.min -= current_flow;
            g_[e].cap_info.max -= current_flow; // TODO these 4
            g_[e].cap_info.avg -= current_flow;
            auto src = boost::source(e, g_);
            auto trg = boost::target(e, g_);
            std::cerr << "Flow remaining on edge " << g_[src].index << "-->" << g_[trg].index << " (" << g_[e].capacity << " of " << g_[e].capacity + current_flow << ")" << std::endl;
            std::cerr << "Corresponding to " << g_[e].capacity/(g_[e].capacity + current_flow) * 100 << "%" << std::endl;
        }
        if (current_flow < min_flow)
            min_flow = current_flow;
    }
    return std::make_pair(contig,std::make_pair(flow_frac, min_flow)); //TODO
}

void UnitigGraph::blockPath(UEdge curr, unsigned int visits)
{
    while (true)
    {
        g_[curr].visits.visits.push_back(visits);
        g_[curr].last_visit = visits;
        auto target = boost::target(curr, g_);
        auto out_edges = boost::out_edges(target, g_);
        float max = 0;
        float max_unvisited = 0;
        UEdge max_e;
        UEdge max_unvisited_e = curr;
        for (auto&& e : out_edges)
        {
            if (g_[e].visits.visits.size() == 0) // an unvisited edge remains
            {
                if (g_[e].cap_info.first > max_unvisited) // is maximal unvisited edge
                {
                    max_unvisited = g_[e].cap_info.first;
                    max_unvisited_e = e;
                }
            }
            else if (g_[e].visits.visits.size() != 0 and g_[e].cap_info.first > max) // maximal visited edge
            {
                max = g_[e].cap_info.first;
                max_e = e;
            }
        }
        if (max != 0 and g_[max_e].last_visit == visits) // the max edge has been visited in the same run
        { // TODO we still might want to continue, if the next edge has a similar coverage as the first had
            return;
        }
        if (max_unvisited != 0) // there still is an unvisited edge
        {
            // check whether we want to continue
            curr = max_unvisited_e;
        }
        else if (max != 0) // all out_edges are visited
        {
            //check whether we want to continue
            curr = max_e;
        }
        else // no out_edges from this edge, check somewhere else
        {
            return;
        }
    }
}

std::pair<UEdge, std::pair<bool, unsigned int>> UnitigGraph::getUnvisitedEdge(const std::vector<UEdge>& sources, unsigned int visits, unsigned int source)
{
    unvisit(); // make sure we havent screwed up
    UEdge curr;
    bool unblocked = false;
    if (visits <= sources.size()) // first check all sources
    {
        unblocked = true;
        curr = sources[visits - 1]; // visits start at 1
    }    
    else // if they have been checked, check whether there are unvisited edges and continue
    { // find the "first one", starting from a source, find the first unvisited edge
        if (sources.size() == 0) // there is no source -> we are complete cycle -> pick random unvisited edge
        {
            for (auto e : boost::edges(g_))
            {
                if (g_[e].last_visit == 0)
                {
                    unblocked = true;
                    curr = e;
                }
            }
        }
        else // at least one source which has been visited before
        {
            auto current_source = sources[source];
            auto target = boost::target(current_source, g_);
            auto out_edges = boost::out_edges(target, g_);
            //breadth first search for the first unvisited edge
            std::queue<UEdge> to_check;
            for (auto e : out_edges)
            {
                to_check.push(e);
            }
            while (!to_check.empty())
            {
                auto next = to_check.front();
                to_check.pop();
                if (g_[next].visited) // to prevent cycling/multiple searches
                    continue;
                else
                    g_[next].visited = true;
                if (g_[next].last_visit == 0)
                {
                    unblocked = true;
                    curr = next;
                    break; //we found an unblocked edge, return this one
                }
                else // edge is visisted, continue search
                {
                    target = boost::target(next, g_);
                    out_edges = boost::out_edges(target, g_);
                    for (auto e : out_edges)
                    {
                        to_check.push(e);
                    }
                }
            }
            if (!unblocked)
            {
                source++;
                if (source < sources.size())
                    return getUnvisitedEdge(sources, visits, source);
                else
                    return std::make_pair(curr, std::make_pair(false, source));
                // we didnt find an unblocked vertex starting from this source, search the next source
            }
        }
    }
    return std::make_pair(curr, std::make_pair(unblocked, source));
}

void UnitigGraph::fixFlow()
{
    unvisit(); // just make sure nothing is screwed up
    std::vector<UEdge> sources;
    for (auto e : boost::edges(g_))
    {
        auto src = boost::source(e, g_);
        auto in_degree = boost::in_degree(src, g_);
        if (in_degree == 0 and std::find(sources.begin(), sources.end(), e) == sources.end()) // only add all source edges if they havent beed added before
        {
            for (auto f : boost::out_edges(src, g_))
            {
                sources.push_back(f); // for every outedge of source start search
            }
        }
    }
    auto edge_compare = [&](UEdge e1, UEdge e2){ //sort by biggest capacity
        return g_[e1].capacity > g_[e2].capacity;
    };
    std::sort(sources.begin(), sources.end(), edge_compare); // so that we search the highest source first
    unsigned int visits = 1;
    unsigned int source = 0;
    while (true)
    {
        std::pair<UEdge, std::pair<bool, unsigned int>> unvisited = getUnvisitedEdge(sources, visits, source);
        UEdge curr = unvisited.first;
        bool unblocked = unvisited.second.first;
        source = unvisited.second.second;
        if (unblocked)
        {
            blockPath(curr, visits++); //marks the first path
        }
        else
            break;
    }
    // we now have the tentative paths, now check how many edges are unique per path
    std::vector<std::vector<UEdge>> unique;
    std::vector<std::vector<UEdge>> total;
    for (unsigned int i = 0; i < visits; i++)
    {
        unique.push_back(std::vector<UEdge>{});
        total.push_back(std::vector<UEdge>{});
    }
    for (auto e : boost::edges(g_))
    {
        for (unsigned int i = 0; i < visits; i++)
        {
            if (g_[e].visits.visits.front() == i + 1)
                unique[i].push_back(e);
            if (std::find(g_[e].visits.visits.begin(), g_[e].visits.visits.end(), i + 1) != g_[e].visits.visits.end())
                total[i].push_back(e);
        }
    }
    std::cerr << "Total edges: " << boost::num_edges(g_) << std::endl;
    for (unsigned int i = 0; i < visits; i++)
    {
        std::cerr << i+1 << " unique/total: " << unique[i].size() << "/" << total[i].size() << std::endl;
        /*if (unique[i].size() < 0.02 * boost::num_edges(g_))
        {
            std::cerr << "Path " << i+1 << " is erroneous" << std::endl;
            std::cerr << "Edges: " << std::endl;
            for (auto e : unique[i])
            {
                std::cerr << g_[e].name << " (Cap: " << g_[e].capacity << ")" << std::endl;
            }
        }*/
    }
}

// calculates the flows and corresponding paths through the graph
void UnitigGraph::assemble(std::string fname)
{
    std::cerr << "Calculating flow..." << std::endl;

    unsigned int i = 0;
    std::cerr << "Contracting simple paths" << std::endl;
    contractPaths();
    /*seed =*/ fixFlow();
    while (true)
    {
        cleanGraph();
        if (boost::num_vertices(g_) == 0)
            break;
        auto seed = getSeed();
        if (g_[seed].capacity <= threshold_)
        {
            break; // all flow has been used, this CC is clear
        }
        std::cerr << "Fixing flow" << std::endl;
        std::cerr << "Finding fattest path..." << std::endl;
        std::vector<UEdge> path = find_fattest_path(seed);
        std::string filename = fname + "Graph" + std::to_string(i) + ".dot";
        std::ofstream outfile (filename);
        printGraph(outfile);
        auto contig = calculate_contigs(path);
        if (contig.first.size() > 150)
        {
            std::cout << ">Contig_" << i++ << "_flow_" << contig.second.first << "_of_" << contig.second.second << std::endl;
            std::cout << contig.first << std::endl;
        }
        else
        {
            std::cerr << "Removed short contig" << std::endl;
            std::cerr << contig.first << std::endl;
            i++;
        }
    }
}

void UnitigGraph::printGraph(std::ostream& os) const
{
    typedef std::map<UVertex, int> IndexMap;
    IndexMap mapIndex;
    boost::associative_property_map<IndexMap> propmapIndex(mapIndex);
    uvertex_iter vi, vi_end;
    int i = 1;
    for (boost::tie(vi,vi_end) = boost::vertices(g_); vi != vi_end; ++vi)
    {
        boost::put(propmapIndex,*vi,i++);
    }

    //boost::write_graphviz(os, g_, boost::make_label_writer(boost::get(&VertexProperties::index,g_)), boost::make_label_writer(boost::get(&EdgeProperties::name,g_)), boost::default_writer(), propmapIndex);
    //boost::write_graphviz(os, g_, boost::make_label_writer(boost::get(&VertexProperties::index,g_)), boost::make_label_writer(boost::get(&EdgeProperties::cap_info,g_)), boost::default_writer(), propmapIndex);
    //boost::write_graphviz(os, g_, boost::make_label_writer(boost::get(&VertexProperties::index,g_)), boost::make_label_writer(boost::get(&EdgeProperties::capacity,g_)), boost::default_writer(), propmapIndex);
    boost::write_graphviz(os, g_, boost::make_label_writer(boost::get(&VertexProperties::index,g_)), boost::make_label_writer(boost::get(&EdgeProperties::visits,g_)), boost::default_writer(), propmapIndex);
    //boost::write_graphviz(os, g_, boost::make_label_writer(boost::get(&VertexProperties::index,g_)), boost::make_label_writer(boost::get(&EdgeProperties::distance2,g_)), boost::default_writer(), propmapIndex);
}

void UnitigGraph::debug()
{
    //markCycles();

	// DEBUG
	auto&& numV = boost::num_vertices(g_);
	auto&& numE = boost::num_edges(g_);
	
	std::cerr << "Unitig graph has " << numV << " vertices and " << numE << " edges, starting cleaning" << std::endl;
}
