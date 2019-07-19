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
UnitigGraph::UnitigGraph() : cc_(1)
{
    std::vector<UVertex> vertices;
    for (unsigned int i = 1; i < 7; i++)
    {
        unsigned int index = 1;
        vertices.push_back(addVertex(&index, std::to_string(i), 1));
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
UnitigGraph::UnitigGraph(deBruijnGraph& dbg, float error_rate) : cc_(1)
{
	std::cerr << "deBruijnGraph has " << dbg.getSize() << " vertices" << std::endl;
	std::cerr << "Building unitig graph from deBruijn graph..." << std::endl;
    
	clock_t t = clock();
	auto&& junc = dbg.getJunctions();
	unsigned int index = 1;
	auto&& out_unbalanced = junc.first;
	auto&& in_unbalanced = junc.second;
    auto&& thresholds = calculate_thresholds(dbg, error_rate);
    thresholds_ = thresholds;
    // starting from the sources, we build the unitig graph
	for (auto& v : out_unbalanced)
	{
		std::string curr = v.get_kmer();
		Vertex* source = dbg.getVertex(curr);
        unsigned int cc = source->cc - 1; // cc start at 1
        float threshold = thresholds[cc];
		// make a guess whether we are relevant already 
		if (!source->is_flagged() and threshold > 0 and source->get_total_out_coverage() > threshold)
		{
			connectUnbalanced(source, &index, curr, dbg, error_rate, threshold); 
			cc_++; // new connected component found
		}
        else
        {
            source->flag();
            //std::cerr << "Skipped small/erroneous dBg" << std::endl;
        }
	}
	for (auto& v : in_unbalanced)
	{
		std::string curr = v.get_kmer();
		Vertex* source = dbg.getVertex(curr);
        unsigned int cc = source->cc - 1;
        float threshold = thresholds[cc];
		if (!source->is_flagged() and threshold > 0 and source->get_total_in_coverage() > threshold)
		{
			connectUnbalanced(source, &index, curr, dbg, error_rate, threshold); 
			cc_++; // new connected component found
		}
        else
        {
            source->flag();
            //std::cerr << "Skipped small/erroneous dBg" << std::endl;
        }
	}
	std::cerr << "Unitig graph successfully build in " << (clock() - t)/1000000. << " seconds." << std::endl;
    std::cerr << "Unitig graph has " << boost::num_vertices(g_) << " vertices" << std::endl;
}

std::vector<float> UnitigGraph::calculate_thresholds(deBruijnGraph& dbg, float error_rate)
{
    std::cerr << "Getting connected components" << std::endl;
    auto t = clock();
    auto dbgs = dbg.split_ccs();
    //delete dbg; // we split it up and can delete the original
    std::cerr << "Getting CCs took " << (clock() - t)/1000000. << " seconds" << std::endl;
    t = clock();
    std::cerr << "Calculating coverage distribution" << std::endl;
    auto&& cov_distr = dbg.coverageDistribution(dbgs);
    std::cerr << "Calculating coverage distribution took " << (clock() - t)/1000000. << " seconds" << std::endl;
    std::vector<float> thresholds;
    for (auto& covs : cov_distr)
    {
        if (covs.empty())
        {
            continue;
        }
        unsigned int members = 0;
        std::vector<float> sorted_coverage;
        sorted_coverage.resize(covs.rbegin()->first + 1, 0.f); // get last (= biggest) element of map
        for (auto&& cov : covs)
        {
            auto pos = cov.first;
            auto val = cov.second;
            members += val;
            sorted_coverage[pos] = val;
        }
        if (members < 500) //less than 500 kmers
        {
            //std::cerr << "^skipped" << std::endl;
            thresholds.push_back(0.f);
            continue;
        }
        //std::cerr << "Coverages: " << std::endl;
        //std::cerr << sorted_coverage << std::endl;
        auto roll = rolling(sorted_coverage, 5); // TODO set window size?
        auto cumm = cummin(roll);
        unsigned int counter = 0;
        unsigned int j = 0;
        for (auto zip : boost::combine(cumm, roll))
        {
            float cummin_val;
            float rolling_val;
            boost::tie(cummin_val, rolling_val) = zip;
            if (cummin_val < rolling_val)
            {
                counter++;
            }
            else
            {
                counter = 0;
            }
            if (counter == 6) //TODO set value (window_size + 1 makes sense)
            {
                std::cerr << "Threshold set to: " << float(j) << std::endl;
                thresholds.push_back(float(j));
                break;
            }
            j++; //position
        }
        if (counter != 6) // no break was encountered (TODO: value)
        {
            //std::cerr << "No signal, threshold set to 1" << std::endl;
            thresholds.push_back(1.f); // no signal found TODO
        }
    }
    return thresholds;
}

// calculates cumulative minimum of in
std::vector<float> UnitigGraph::cummin(std::vector<float>& in)
{
    std::vector<float> cummin;
    cummin.reserve(in.size());
    float curr = in.front();
    for (auto& val : in)
    {
        if (std::isnan(curr) or (val < curr and !std::isnan(val)))
        {
            curr = val;
        }
        cummin.push_back(curr);
    }
    return cummin;
}

// calculates rolling average over window size "len"
std::vector<float> UnitigGraph::rolling(std::vector<float>& in, unsigned int len)
{
    std::vector<float> rolling;
    rolling.resize(in.size());
    float total = 0.0;
    for (unsigned int i = 0; i < in.size() ; i++)
    {
        total += float(in.at(i));
        if (i >= len)
        {
            total -= in.at(i - len);
        }
        if (i >= (len - 1))
        {
            rolling.at(i) = total/float(len);
        }
        else
        {
            rolling.at(i) = std::nanf(""); //add nan values for first len - 1 elements
        }
    }
    return rolling;
}

// adds a vertex to the unitig graph: adds it to the boost graph, as well as to the mapping from index to vertex
UVertex UnitigGraph::addVertex(unsigned int* index, std::string name, unsigned int ccc)
{
	UVertex uv = boost::add_vertex(g_);
	(*index)++;
    // set vertex properties
    g_[uv].name = name;
	g_[uv].index = *index;
    g_[uv].scc = 1;
    g_[uv].tarjan_index = 0; // needs to be 0 to find out whether it has been set
    g_[uv].onStack = false;
    g_[uv].cc = ccc - 1; // set cc to current CC (starting with 0 here)
    g_[uv].visiting_time = 0;
	
    auto&& ins = std::make_pair(*index,uv);
	graph_.insert(ins);
	return uv;
}

// function for connecting a given source/sink vertex to all its unbalanced successors/predecessors
void UnitigGraph::connectUnbalanced(Vertex* source, unsigned int* index, std::string curr, deBruijnGraph& dbg, float error_rate, float threshold)
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
            if (junction->get_out_coverage(c) > error_rate * junction->get_total_out_coverage()/* and junction->get_out_coverage(c) > threshold_*/)
			{
				to_search = true; break;
			}
		}
		if (!to_search)
		{
			for (const char& c : pred)
			{
				if (junction->get_in_coverage(c) > error_rate * junction->get_total_in_coverage()/* and junction->get_in_coverage(c) > threshold_*/)
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
			uv = addVertex(index, seq, junction->cc);
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
		auto&& following = addNeighbours(seq, succ, pred, dbg, index, uv, threshold); // finding the next unbalanced vertices
		for (auto v : following)
		{	
			// if no sequence is returned, no vertex was added, so we do not need to continue on this vertex
			if (v.second != "")
				todo.push(v);
		}
	} while (!todo.empty());
}

// iterating over all neighbours of current node, build the different sequences to the next unbalanced node
std::vector<std::pair<Vertex*,std::string> > UnitigGraph::addNeighbours(std::string& curr, const std::vector<char>& succ, const std::vector<char>& pred, deBruijnGraph& dbg, unsigned int* index, UVertex& uv, float threshold)
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
                float curr_start = 0.f;
                if (currV->get_total_out_coverage() != 0)
                    curr_start = currV->get_read_starts() * currV->get_out_coverage(n)/currV->get_total_out_coverage();
                next = curr.substr(1) + n;
                Vertex* nextV = dbg.getVertex(next);
                Sequence s = *dbg.getSequence(next);
                sequence += n;
                if (!nextV->is_flagged())
                    following.push_back(buildEdge(uv, nextV, next, sequence, index, coverage, pcov, dbg, curr_start, threshold));
            }
            // if we are a reverse complement, we actually want to add the path in reverse order
            else
            {
                float curr_end = 0.f;
                if (currV->get_total_in_coverage() != 0)
                    curr_end = currV->get_read_ends() * currV->get_in_coverage(deBruijnGraph::complement(n))/currV->get_total_in_coverage();
                next = deBruijnGraph::complement(n) + curr.substr(0,curr.length() - 1);
                Vertex* nextV = dbg.getVertex(next);
                Sequence s = *dbg.getSequence(next);
                sequence += curr.back(); // the predecessor points to the current vertex with the last char of curr (by definition)
                if (!nextV->is_flagged())
                    following.push_back(buildEdgeReverse(uv, nextV, next, sequence, index, coverage, pcov, dbg, curr_end, threshold));
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
                float curr_end = 0.f;
                if (currV->get_total_in_coverage() != 0)
                    curr_end = currV->get_read_ends() * currV->get_in_coverage(n)/currV->get_total_in_coverage();
                prev = n + curr.substr(0,curr.length() - 1);
                Vertex* nextV = dbg.getVertex(prev);
                Sequence s = *dbg.getSequence(prev);
                sequence += curr.back(); // the predecessor points to the current vertex with the last char of curr
                if (!nextV->is_flagged())
                    following.push_back(buildEdgeReverse(uv, nextV, prev, sequence, index, coverage, pcov, dbg, curr_end, threshold));
            }
            else
            {
                float curr_start = 0.f;
                if (currV->get_total_out_coverage() != 0)
                    curr_start = currV->get_read_starts() * currV->get_out_coverage(deBruijnGraph::complement(n))/currV->get_total_out_coverage();
                prev = curr.substr(1) + deBruijnGraph::complement(n);
                Vertex* nextV = dbg.getVertex(prev);
                Sequence s = *dbg.getSequence(prev);
                sequence += deBruijnGraph::complement(n);
                if (!nextV->is_flagged())
                    following.push_back(buildEdge(uv, nextV, prev, sequence, index, coverage, pcov, dbg, curr_start, threshold));
            }
        }
    }
	currV->flag(); // this vertex is done
	return following;
}

// go back through the graph until the next unbalanced node is found and add an ("reversed") edge
std::pair<Vertex*,std::string> UnitigGraph::buildEdgeReverse(UVertex trg, Vertex* nextV, std::string prev, std::string& sequence, unsigned int* index, float coverage, float pcov, deBruijnGraph& dbg, float curr_end, float threshold)
{
	char lastchar = g_[trg].name.back(); //char with which we are pointing to trg
    float total_out = nextV->get_total_out_coverage();
    if (total_out == 0)
    {
        unsigned int comp = 1;
        total_out = std::max(nextV->get_out_coverage(lastchar), comp);
    }
     // total_out is >= 1
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
	}
	avg /= float(length); // average coverage over the path
	if (!nextV->is_visited() and (avg >= threshold or sequence.length() > 500)) // TODO if coverage is low but the (unique) sequence is long, still add
	{// if the next vertex has been visited it already is part of the unitiggraph, otherwise add it
		nextV->visit();
		addVertex(index, prev, nextV->cc); // the vertex is new and found to be relevant
		nextV->index = *index;
	}
	else if (!nextV->is_visited() or (avg < threshold and sequence.length() <= 500) or nextV->index == 0)
	{
		return std::make_pair(nextV,""); // path has too low coverage
	}
	UVertex src = graph_[nextV->index];
	auto e = boost::edge(src,trg,g_);
	// if edge has been added or the immediate neighbour is an unbalanced vertex, do not add edge
	if ((!e.second or (e.second and g_[e.first].name != sequence)))
	{
        //set new edge's information
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
        g_[e.first].residual_capacity = avg;
        g_[e.first].cap_info.avg = avg;
        g_[e.first].cap_info.max = max;
        g_[e.first].cap_info.min = min;
        g_[e.first].cap_info.first = first;
        g_[e.first].cap_info.last = last;
        g_[e.first].cap_info.length = g_[e.first].name.length();
        if (g_[e.first].name.length() > 0 and avg > 0)
        {
            g_[e.first].cap_info.starting = starts_with/(g_[e.first].cap_info.length * avg);
            g_[e.first].cap_info.ending = ends_with/(g_[e.first].cap_info.length * avg);
        }
        else // shouldn't end up here
        {
            g_[e.first].cap_info.starting = starts_with;
            g_[e.first].cap_info.ending = ends_with;
        }
        g_[e.first].visited = false;
        g_[e.first].first_vertex = false;
	}
	return std::make_pair(nextV,prev);
}

// same function like the reverse one, but going forward and finding successors
std::pair<Vertex*,std::string> UnitigGraph::buildEdge(UVertex src, Vertex* nextV, std::string next, std::string& sequence, unsigned int* index, float coverage, float pcov, deBruijnGraph& dbg, float curr_start, float threshold)
{
	// with a little effort this can be moved inside the while loop for efficiency reasons
    float first_char = g_[src].name.front();
    float starts_with = curr_start; // number of reads starts within this edge, normalized by total flow through this edge
    float total_in = nextV->get_total_in_coverage();
    if (total_in == 0)
    {
        unsigned int comp = 1;
        total_in = std::max(nextV->get_in_coverage(first_char), comp);
    } // total_in >= 1
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
	if (!nextV->is_visited() and (avg >= threshold or sequence.length() > 500)) //TODO arbitrary value
	{
		nextV->visit();
		addVertex(index, next, nextV->cc);
		nextV->index = *index;
	}
	else if (!nextV->is_visited() or (avg < threshold and sequence.length() <= 500) or nextV->index == 0)
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
        g_[e.first].residual_capacity = avg;
        g_[e.first].cap_info.avg = avg;
        g_[e.first].cap_info.max = max;
        g_[e.first].cap_info.min = min;
        g_[e.first].cap_info.first = first;
        g_[e.first].cap_info.last = last;
        g_[e.first].cap_info.length = g_[e.first].name.length();
        if (g_[e.first].name.length() > 0 and avg > 0)
        {
            g_[e.first].cap_info.starting = starts_with/(g_[e.first].cap_info.length * avg);
            g_[e.first].cap_info.ending = ends_with/(g_[e.first].cap_info.length * avg);
        }
        else // shouldn't end up here
        {
            g_[e.first].cap_info.starting = starts_with;
            g_[e.first].cap_info.ending = ends_with;
        }
        g_[e.first].visited = false;
        g_[e.first].first_vertex = false;
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
		unsigned int outdegree = boost::out_degree(*vi, g_);

		// if in and outdegree is 1, we are on a simple path and can contract again
		if (outdegree == 1 and indegree == 1)
		{
			auto&& ie = boost::in_edges(*vi,g_);
			auto&& oe = boost::out_edges(*vi,g_);
			auto&& new_source = boost::source(*ie.first,g_);
			auto&& new_target = boost::target(*oe.first,g_);
            if (new_source == new_target)
            {
                continue; // do not contract to single vertex (which might get deleted)
            }
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
            if (seq.length() > 0)
			    capacity /= seq.length(); // currently using the average coverage on the contracted path

			auto&& new_e = boost::add_edge(new_source,new_target,g_);
            g_[new_e.first].last_visit = 0;
            g_[new_e.first].name = seq;
            g_[new_e.first].capacity = capacity;
            g_[new_e.first].residual_capacity = capacity;
            g_[new_e.first].cap_info.avg = capacity;
            g_[new_e.first].cap_info.max = max;
            g_[new_e.first].cap_info.min = min;
            g_[new_e.first].cap_info.first = first;
            g_[new_e.first].cap_info.last = last;
            g_[new_e.first].cap_info.length = g_[new_e.first].name.length();
            g_[new_e.first].visited = false;
            g_[new_e.first].first_vertex = false;
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


bool UnitigGraph::hasRelevance()
{
    unsigned int length = 0;
    for (auto&& e : boost::edges(g_))
    {
        length += g_[e].name.size();
    }
    if (length <= 500)
    {
        std::cerr << "Graph undercutting threshold of 500 characters (" << length << ")" << std::endl;
    }
    return length > 500;
}

// the graph might contain some unconnected vertices, clean up
void UnitigGraph::cleanGraph()
{
    for (auto e : boost::edges(g_))
    {
        g_[e].last_visit = 0; // reusing for number of allowed paths
    }
	removeEmpty();
	removeStableSets();
}

// Tests whether two percentages "belong together"
bool UnitigGraph::test_hypothesis(float to_test_num, float to_test_denom, float h0, float threshold)
{
    float diff = std::abs(to_test_num - h0 * to_test_denom); // the absolute difference between num and expected num
    //bool less = to_test_num < to_test_denom * h0;
    //float perc = to_test_num/to_test_denom; //currently unused
    return (/*less or*/ diff < threshold); // they differ by less than threshold 
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

// run dijsktra with fatness as optimality criterion, marks the graph with the distances from seed
void UnitigGraph::dijkstra(UEdge seed, bool init)
{
    auto edge_compare = [&](UEdge e1, UEdge e2){ //sort by biggest fatness
        return g_[e1].fatness < g_[e2].fatness;
    };
    std::vector<UEdge> q;
    if (init)
    {
        for (auto e : boost::edges(g_)) //initialise distances and fatness
        {
            g_[e].fatness = 0;
            g_[e].distance = std::numeric_limits<unsigned int>::max();
            if (e == seed) //
            {
                g_[e].distance = 0;
                g_[e].fatness = std::numeric_limits<float>::max();
                // such that source plays a lesser role
            }
        }
    }
    if (!init) // if fatness is not initialised for all vertices, initialize for all following
    {
        q.push_back(seed);
        while (!q.empty()) // classic dijsktra routine (cancelling when cycle found)
        {
            auto curr = q.back();
            q.pop_back();
            g_[curr].visited = true;
            g_[curr].fatness = 0;
            auto target = boost::target(curr, g_);
            for (auto oe : boost::out_edges(target, g_))
            {
                if (!g_[oe].visited)
                {
                    q.push_back(oe);
                }
            }
        }
    }
    g_[seed].fatness = std::numeric_limits<float>::max();
    if (g_[seed].distance == std::numeric_limits<float>::max())
    {
        g_[seed].distance = 0;
    }
    for (auto e : boost::edges(g_)) //unvisit for next run
    {
        g_[e].visited = false;
    }
    q.push_back(seed);
    while (!q.empty()) // classic dijsktra routine (cancelling when cycle found)
    {
        auto curr = q.back();
        q.pop_back();
        g_[curr].visited = true;
        auto target = boost::target(curr, g_);
        for (auto oe : boost::out_edges(target, g_))
        {
            float fat = g_[oe].fatness;
            if (fat < std::min(g_[curr].fatness, (init ? g_[oe].residual_capacity : g_[oe].capacity)))
            {
                g_[oe].fatness = std::min(g_[curr].fatness, (init ? g_[oe].residual_capacity : g_[oe].capacity));
                g_[oe].prev = curr;
                g_[oe].distance = g_[curr].distance + g_[oe].name.size();
            }
            if (!g_[oe].visited)
            {
                q.push_back(oe);
            }
            else
            {
                g_[oe].first_vertex = true; // this is the first vertex of a cycle
            }
        }
        std::sort(q.begin(), q.end(), edge_compare);
    }
    g_[seed].fatness = (init ? g_[seed].residual_capacity : g_[seed].capacity);
    for (auto e : boost::edges(g_)) //unvisit for next run
    {
        g_[e].visited = false;
    }
}

std::pair<UEdge, float> UnitigGraph::get_target(UEdge seed, bool lenient)
{
    float max_dist = 0;
    float running_fatness = 0;
    float running_distance = 0;
    UEdge last;
    auto visits = g_[seed].visits;
    //bool first_vertex = false;
    unvisit();
    std::vector<UEdge> q = {seed};
    while (!q.empty()) // classic dijsktra routine (cancelling when cycle found)
    {
        auto curr = q.back();
        q.pop_back();
        
        bool same_visit = false;
        
        for (auto v : g_[curr].visits)
        {
            same_visit = (std::find(visits.begin(), visits.end(), v) != visits.end());
            if (same_visit)
                break;
        }

        g_[curr].visited = true;
        unsigned int num_visits = g_[curr].visits.size();
        auto target = boost::target(curr, g_);
        /*if (same_visit and num_visits == 1 and g_[curr].distance < std::numeric_limits<unsigned int>::max() and (running_distance == 0 or (g_[curr].fatness > 0 and running_fatness/g_[curr].fatness < g_[curr].distance/running_distance)))
        {
            running_fatness = g_[curr].fatness;
            running_distance = g_[curr].distance;
            last = curr;
            max_dist = g_[curr].distance;
        }*/
        if (g_[curr].distance > running_distance and g_[curr].distance != std::numeric_limits<unsigned int>::max() and same_visit)
        {
            running_distance = g_[curr].distance;
            last = curr;
            max_dist = g_[curr].distance;
        }
        for (auto oe : boost::out_edges(target, g_))
        {
            if (!g_[oe].visited)
            {
                q.push_back(oe);
            }
        }
    }
    if (max_dist == 0)
    {
        last = seed; // so we don't return nothing
    }
    unvisit();
    return std::make_pair(last, max_dist);
}

std::vector<UEdge> UnitigGraph::fixFlow(UEdge seed)
{
    unvisit();
    dijkstra(seed, false);
    auto path = find_fattest_path(seed);
    float flow = 0.;
    unsigned int edges = 0;
    bool corrected = true;
    for (auto& e : path)
    {
        if (g_[e].visits.size() == 1)
        {
            flow += g_[e].capacity;
            edges++;
        }
    }
    if (edges > 0) // if unique  edges are present calculate average flow over them
    {
        flow /= edges;
    }
    else // else take flow of first edge of path (TODO?)
    {
        flow = g_[path.front()].capacity;
    }
    unsigned int pos = 0;
    while (corrected)
    {
        bool eq = false;
        unsigned int i = 0;
        for (auto& e : path)
        {
            auto trg = boost::target(e, g_);
            float val = -1;
            for (auto&& oe : boost::out_edges(trg, g_))
            {
                if (g_[oe].fatness == val)
                {
                    eq = true;
                }
                else
                {
                    val = g_[oe].fatness;
                }
            }
            if (eq and i > pos)
            {
                unvisit();
                pos = i;
                dijkstra(e, false);
                break;
            }
            i++;
        }
        corrected = (eq and i == pos);
    }
    unvisit();
    std::vector<UEdge> return_path = find_fattest_path(seed);
    return return_path;
}

// Calculates the fattest path through the graph, marks vertices as being visited
std::vector<UEdge> UnitigGraph::find_fattest_path(UEdge seed)
{
    auto source = boost::source(seed, g_);
    //auto target = boost::target(seed, g_);
    //std::cerr << "Source: " << g_[source].index << " -> " << g_[target].index << ": " << g_[seed].capacity << std::endl;
    auto trg = get_target(seed, false);
    auto last = trg.first;
    //auto source2 = boost::source(last, g_);
    //auto target2 = boost::target(last, g_);
    //std::cerr << "Target: " << g_[source2].index << " -> " << g_[target2].index << ": " << g_[last].capacity << std::endl;
    float max_dist = trg.second;
    if (max_dist == 0) // path is only one edge, no longest path
    {
        return std::vector<UEdge>{seed};
    }
    auto curr = last;
    std::deque<UEdge> path = {curr};
    g_[curr].last_visit++;
    while (curr != seed and g_[curr].distance > 0 and g_[curr].distance < std::numeric_limits<unsigned int>::max() and !g_[curr].visited) // this means the distance has been set, i.e. the vertex has been reached
    {
        g_[curr].visited = true;
        curr = g_[curr].prev;
        g_[curr].last_visit++;
        path.push_front(curr);
    }
    source = boost::source(path.front(), g_);
    unsigned int i = 0;
    unsigned int j = 1;
    float seq_length = 0;
    std::vector<UEdge> ret;
    
    //std::cerr << g_[source].index; 
    for (auto e : path)
    {
        auto trg = boost::target(e, g_);
        //std::cerr << " -> " << g_[trg].index << "(" << g_[e].capacity << " " << g_[e].fatness << ") ";
        seq_length += g_[e].name.size();
        auto ct = 0;
        auto avg = 0.f;
        for (unsigned int k = i; k < j; k++)
        {
            unsigned int visits = g_[path[k]].visits.size();
            if (visits < 2)
            {
                avg += g_[path[k]].capacity;
                ct++;
            }
        }
        if (ct == 0 or avg == 0)
        {
            ret.push_back(e);
            continue;
        }
        avg /= ct;
        if (seq_length > 500)
        {
            seq_length -= g_[path[i]].name.size();
            i++;
        }
        j++;
        ret.push_back(e);
    }
    //std::cerr << std::endl;
    return ret;
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
        g_[e].first_vertex = false;
    }
}

void UnitigGraph::reduce_flow(std::vector<UEdge>& path, float flow, std::vector<float>& flows, std::set<unsigned int>& unique_paths)
{
    float removed_coverage = g_[path.front()].capacity;
    for (auto e : path)
    {
        // first find out which path we are on (we delete this because it has been used then)
        auto to_remove = g_[e].visits.begin();
        for (auto p : unique_paths)
        {
            to_remove = std::find(g_[e].visits.begin(), g_[e].visits.end(), p);
            if (to_remove != g_[e].visits.end())
            {
                break; // TODO
            }
        }
        // TODO fix multiple "unique" paths
        if (to_remove != g_[e].visits.end())
        {
            g_[e].visits.erase(to_remove);
        }
        else
        {
            if (!g_[e].visits.empty())
            {
                g_[e].visits.erase(g_[e].visits.begin());
            }
        }
        // now check whether this was the last path or there are remaining paths to reduce capacity
        float val = g_[e].capacity;
        float val_last = g_[e].cap_info.last;
        float val_first = g_[e].cap_info.first;
        auto cc = g_[boost::source(e, g_)].cc;
        auto threshold = thresholds_[cc];
        
        if (g_[e].visits.empty())
        {
            g_[e].capacity = 0;
            g_[e].cap_info.first = 0;
            g_[e].cap_info.last = 0;
            g_[e].cap_info.avg = 0;
            removed_coverage = val;
        }
        else
        {
            bool not_decreasing = g_[e].cap_info.first <= 1.1 * g_[e].cap_info.last or std::abs(g_[e].cap_info.first - g_[e].cap_info.last) < threshold; //test_hypothesis(g_[e].cap_info.first, g_[e].cap_info.last, 1.2, threshold_); //not decreasing (first/last >= 1.2)
            bool not_increasing = g_[e].cap_info.last <= 1.1 * g_[e].cap_info.first or std::abs(g_[e].cap_info.first - g_[e].cap_info.last) < threshold; //test_hypothesis(g_[e].cap_info.last, g_[e].cap_info.first, 1.2, threshold_); // not increasing (last/first >= 1.2)
            //cannot both be false, but both be true (if close to 1.2)
            if (val > threshold and not_decreasing and not_increasing)
            {
                g_[e].capacity = std::max(threshold, g_[e].capacity - removed_coverage); // there might be paths, so leave a small amount
                g_[e].cap_info.first = std::max(threshold, g_[e].cap_info.first - removed_coverage);
                g_[e].cap_info.last = std::max(threshold, g_[e].cap_info.last - removed_coverage);
                g_[e].cap_info.avg = g_[e].capacity;
                removed_coverage = val - g_[e].capacity;
            }
            else if (val > threshold and not_decreasing and !not_increasing)
            {
                g_[e].capacity = std::max(threshold, g_[e].cap_info.first - removed_coverage);
                g_[e].cap_info.first = g_[e].capacity;
                g_[e].cap_info.last = std::max(threshold, g_[e].cap_info.last - (val - g_[e].capacity));
                g_[e].cap_info.avg = g_[e].capacity;
                removed_coverage = val_last - g_[e].capacity;
            }
            else if (val > threshold and !not_decreasing and not_increasing)
            {
                g_[e].capacity = std::max(threshold, g_[e].cap_info.last - removed_coverage);
                g_[e].cap_info.last = g_[e].capacity;
                g_[e].cap_info.first = std::max(threshold, g_[e].cap_info.first - (val - g_[e].capacity));
                g_[e].cap_info.avg = g_[e].capacity;
                removed_coverage = val_first - g_[e].capacity;
            }
            else
            {
                g_[e].cap_info.first = 0;
                g_[e].cap_info.last = 0;
                g_[e].capacity = 0;
                g_[e].cap_info.avg = 0;
                removed_coverage = val;
            }
        }
    }
}

// Given all the chosen edges and their coverage fraction, builds the contigs and reduces flow accordingly
std::pair<std::string, float> UnitigGraph::calculate_contigs(std::vector<UEdge>& path, std::vector<float>& flows)
{
    float flow = 0.;
    float max_flow = 0.;
    unsigned int i = 0;
    UEdge curr = path.front();
    UVertex source = boost::source(curr, g_);
    std::set<unsigned int> paths;
    std::string contig = g_[source].name;
    for (auto e : path)
    {
        contig += g_[e].name;
        if (g_[e].visits.size() == 1)
        {
            if (g_[e].capacity > max_flow)
            {
                max_flow = g_[e].capacity;
            }
            paths.insert(g_[e].visits[0]);
            flow += g_[e].capacity;
            i++;
        }
    }
    if (i > 0)
    {
        flow /= i; // average flow over unqiue edges of path
    }
    else // no unique edges on path, remove duplicate path as far as possible and calculate new flow based on removed path
    {
        flow = 0;
        std::vector<unsigned int> to_remove = g_[path.front()].visits; // choose one of these visits to be removed;
        for (auto e : path)
        {
            if (g_[e].visits.size() > 1)
            {
                bool removed = false;
                for (auto v : to_remove)
                {
                    auto visit = std::find(g_[e].visits.begin(), g_[e].visits.end(), v);
                    if (visit != g_[e].visits.end())
                    {
                        g_[e].visits.erase(visit);
                        removed = true;
                        flow += flows[*visit - 1];
                        break;
                    }
                }
                if (!removed)
                {
                    flow += flows[g_[e].visits.front() - 1];
                    g_[e].visits.erase(g_[e].visits.begin());
                }
            }
        }
        std::cerr << "Removed duplicated contig " << std::endl;
        std::cerr << contig << std::endl;
        contig = "";
        if (path.size() > 0)
            flow /= path.size();
    }
    reduce_flow(path, flow, flows, paths);
    return std::make_pair(contig, flow);
}

std::vector<UEdge> UnitigGraph::blockPath(UEdge curr, unsigned int visits)
{
    std::vector<UEdge> blockedPath;
    while (true)
    {
        g_[curr].visits.push_back(visits);
        g_[curr].last_visit = visits;
        blockedPath.push_back(curr);
        auto target = boost::target(curr, g_);
        auto out_edges = boost::out_edges(target, g_);
        float max = -1;
        float max_unvisited = -1;
        UEdge max_e;
        UEdge max_e_unvisited;
        for (auto&& e : out_edges)
        {
            if (g_[e].visits.size() == 0)
            {
                if (g_[e].residual_capacity > max_unvisited)
                {
                    max_unvisited = g_[e].residual_capacity;
                    max_e_unvisited = e;
                }
            }
            else if (g_[e].residual_capacity > max) // maximal visited edge
            {
                max = g_[e].residual_capacity;
                max_e = e;
            }
        }
        if (max != -1 and g_[max_e].last_visit == visits) // the max edge has been visited in the same run
        { // TODO we still might want to continue, if the next edge has a similar coverage as the first had
            return blockedPath;
        }
        if (max_unvisited != -1)
        {
            curr = max_e_unvisited;
        }
        else if (max != -1) // all out_edges are visited
        {
            //check whether we want to continue
            curr = max_e;
        }
        else // no out_edges from this edge, check somewhere else
        {
            return blockedPath;
        }
    }
}

UEdge UnitigGraph::get_next_source() /// just returns the highest capacity edge (TODO?)
{
    auto edge_compare = [&](UEdge e1, UEdge e2){ //sort by biggest capacity
        return g_[e1].capacity > g_[e2].capacity;
    };
    UEdge source;

    auto sources = get_sources();
    if (sources.size() > 0) // if there are sources, take highest possible source
    {
        std::sort(sources.begin(), sources.end(), edge_compare);
        source = sources.front();
    }
    else // else take highest unvisited edge
    {
        float max = 0;
        for (auto e : boost::edges(g_))
        {
            if (g_[e].capacity > max)
            {
                max = g_[e].capacity;
                source = e;
            }
        }
    }
    return source;
}

std::pair<UEdge, bool> UnitigGraph::checkUnvisitedEdges(UEdge current_source)
{
    bool unblocked = false;
    UEdge curr = current_source;
    auto source = boost::source(current_source, g_);
    auto out_edges = boost::out_edges(source, g_);
    //breadth first search for the first unvisited edge
    std::queue<UEdge> to_check;
    to_check.push(curr);
    float capacity = 0.f;
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
        if (g_[next].last_visit == 0) // next edge has not been visited
        {
            if (g_[next].residual_capacity > capacity) // and has higher capacity than other close-by edges
            {
                unblocked = true;
                curr = next;
                capacity = g_[next].residual_capacity; // choose as next
            }
        }
        else // edge is visisted, continue search
        {
            source = boost::target(next, g_);
            out_edges = boost::out_edges(source, g_);
            for (auto e : out_edges)
            {
                if (!unblocked or g_[e].last_visit == 0) // if we havent found an edge or next edge is unvisited: continue search
                    to_check.push(e);
            }
        }
    }
    return std::make_pair(curr, unblocked);
}

std::pair<UEdge, bool> UnitigGraph::getUnvisitedEdge(const std::vector<UEdge>& sources, unsigned int visits)
{
    UEdge curr;
    bool unblocked = false;
    if (sources.size() == 0) // there is no source -> we are complete cycle -> pick highest capacity edge
    {
        float max = -1;
        for (auto e : boost::edges(g_))
        {
            if (g_[e].last_visit == 0 and g_[e].residual_capacity > max)
            {
                unblocked = true;
                curr = e;
                max = g_[e].residual_capacity;
            }
        }
    }
    else // at least one source which has been visited before
    {
        curr = sources.back();
        unblocked = (g_[curr].last_visit == 0);
        for (auto e : sources) // check all sources for unchecked edges
        {
            auto nextUnvisited = checkUnvisitedEdges(e);
            auto potential_source = nextUnvisited.first;
            if (nextUnvisited.second and g_[potential_source].residual_capacity > g_[curr].residual_capacity)
            { // and check the highest capacity one
                unblocked = true;
                curr = potential_source;
            }
        }
    }
    return std::make_pair(curr, unblocked);
}

std::vector<UEdge> UnitigGraph::get_sources()
{
    std::vector<UEdge> sources;
    for (auto e : boost::edges(g_))
    {
        auto src = boost::source(e, g_);
        auto target = boost::target(e, g_);
        auto in_degree = boost::in_degree(src, g_);
        if (in_degree == 0 and std::find(sources.begin(), sources.end(), e) == sources.end()) // only add all source edges if they havent beed added before
        {
            for (auto f : boost::out_edges(src, g_))
            {
                sources.push_back(f); // for every outedge of source start search
            }
        }
        else if (in_degree == 1 and src == target) // also add sources which have a self-loop
        {
            sources.push_back(e);
        }
    }
    return sources;
}

std::vector<float> UnitigGraph::find_paths()
{
    std::vector<UEdge> sources = get_sources(); //get sources of the graph (indegree = 0)
    auto edge_compare = [&](UEdge e1, UEdge e2){ //sort by biggest capacity
        return g_[e1].capacity > g_[e2].capacity;
    };
    std::sort(sources.begin(), sources.end(), edge_compare); // so that we search the highest source first
    unsigned int used_sources = 1;
    unsigned int visits = 1;
    std::vector<std::vector<UEdge>> unique;
    std::vector<UEdge> started_from;
    std::vector<float> unique_paths;
    while (true)
    {
        unvisit();
        std::pair<UEdge, bool> unvisited = getUnvisitedEdge(sources, used_sources);
        UEdge curr = unvisited.first;
        bool unblocked = unvisited.second;
        std::vector<UEdge> blockedPath;
        if (unblocked)
        {
            unique.push_back(std::vector<UEdge>{});
            dijkstra(curr, true);
            started_from.push_back(curr); // add the edge from which we started
            blockedPath = blockPath(curr, visits); //marks the first path
        }
        else
            break;
        float avg = 0;
        unsigned int length = 0;
        for (auto e : boost::edges(g_))
        {
            if (!g_[e].visits.empty() and g_[e].visits.front() == visits)
            {
                unique[visits - 1].push_back(e);
                avg += g_[e].residual_capacity; 
                length += g_[e].name.size();
            }
        }
        auto size = unique[visits - 1].size();
        if (size > 0)
            avg /= size;
        // we now have the tentative paths, now check how many edges are unique per path
        unique_paths.push_back(remove_non_unique_paths(unique, blockedPath, length, visits - 1));
        auto source = boost::source(curr, g_);
        auto target = boost::target(curr, g_);
        auto indegree = boost::in_degree(source, g_);
        if (indegree == 0 or (indegree == 1 and source == target))
        {
            used_sources++;
        }
        visits++;
    }
    return unique_paths;
}

float UnitigGraph::remove_non_unique_paths(std::vector<std::vector<UEdge>>& unique, std::vector<UEdge>& blockedPath, unsigned int length, unsigned int visits)
{
    auto edge_compare = [&](UEdge e1, UEdge e2){ //sort by biggest capacity
        return g_[e1].capacity > g_[e2].capacity;
    };
    float median = 0.f;
    auto size = unique[visits].size();
    if (size < 0.02 * boost::num_edges(g_) and size < 15 and length < 500) //TODO parameters
    {
        for (auto e : blockedPath)
        {
            if (g_[e].visits.size() > 1) // so we dont not visit some edges
            {
                auto pos = std::find(g_[e].visits.begin(), g_[e].visits.end(), visits + 1);
                if (pos != g_[e].visits.end())
                {
                    g_[e].visits.erase(pos);
                }
            }
            else
            {
                g_[e].residual_capacity = 0;
            }
        }
    }
    else
    {
        std::sort(unique[visits].begin(), unique[visits].end(), edge_compare);
        if (visits == 0)
            median = g_[unique[visits][(3*size)/4]].capacity; //first quartile so we get a "unique" edge more likely
        else
            median = g_[unique[visits][size/2]].capacity; //for the second run the median is fine
        for (auto e : unique[visits])
        {
            auto cc = g_[boost::source(e, g_)].cc;
            g_[e].residual_capacity = std::max(thresholds_[cc], g_[e].residual_capacity - median);
        }
    }
    return median;
}

// calculates the flows and corresponding paths through the graph
void UnitigGraph::assemble(std::string fname)
{
    unsigned int i = 0;
    std::cerr << "Cleaning graph" << std::endl;
    cleanGraph();
    contractPaths();
    std::cerr << boost::num_vertices(g_) << " vertices remaining" << std::endl;
    std::cerr << "Calculating paths" << std::endl;
    auto all_paths = find_paths();
    /*DEBUG*/
    //std::string filename = fname + "Graph.dot";
    //std::ofstream outfile (filename);
    //printGraph(outfile);
    /*DEBUG*/
    UEdge seed;
    for (auto e : boost::edges(g_))
    {
        g_[e].last_visit = 0; // reset last visit
    }
    while (hasRelevance())
    {
        unvisit(); // to track visits and residual capacity
        auto seed = get_next_source();
        std::cerr << "Fixing flow and finding fattest path" << std::endl;
        std::vector<UEdge> path = fixFlow(seed);
        std::string filename = fname + "Graph" + std::to_string(i) + ".dot";
        std::ofstream outfile (filename);
        printGraph(outfile);
        std::cerr << "Calculating contig " << i << "..." << std::endl;
        auto contig = calculate_contigs(path, all_paths);
        if (contig.first.size() > 150)
        {
            std::cout << ">Contig_" << i << "_flow_" << contig.second << std::endl;
            std::cout << contig.first << std::endl;
        }
        else
        {
            std::cerr << "Removed short contig" << std::endl;
            std::cerr << contig.first << std::endl;
        }
        i++;
        std::cerr << "Cleaning graph again..." << std::endl;
        cleanPath(path);
    }
    std::cerr << "Assembly complete" << std::endl;
}

void UnitigGraph::removeEmpty()
{
    std::vector<UEdge> toDelete;
    for (auto&& e : boost::edges(g_))
    {
        auto&& src = boost::source(e, g_); //src and target always have the same cc
        if (g_[e].capacity == 0 or g_[e].capacity < thresholds_[g_[src].cc]) //to make sure empty edges are definitely deleted
        {
            toDelete.push_back(e);
        }
    }
    for (auto&& e : toDelete)
    {
        boost::remove_edge(e, g_);
    }
}

void UnitigGraph::cleanPath(std::vector<UEdge>& path)
{
    std::set<UEdge> toDelete;
    for (auto&& e : path)
    {
        auto&& source = boost::source(e, g_);
        if (g_[e].capacity == 0 or g_[e].capacity < thresholds_[g_[source].cc])
        {
            toDelete.insert(e);
        }
    }
    for (auto&& e : toDelete)
    {
        auto source = boost::source(e, g_);
        boost::remove_edge(e, g_);
		unsigned int indegree = boost::in_degree(source, g_);
		unsigned int outdegree = boost::out_degree(source,g_);
		if (outdegree == 0 and indegree == 0)
		{
			boost::remove_vertex(source,g_);
		}
    }
}

void UnitigGraph::printGraph(std::ostream& os)
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
    for (auto e : boost::edges(g_))
    {
        g_[e].v.visits = g_[e].visits;
    }
    //boost::write_graphviz(os, g_, boost::make_label_writer(boost::get(&VertexProperties::index,g_)), boost::make_label_writer(boost::get(&EdgeProperties::name,g_)), boost::default_writer(), propmapIndex);
    boost::write_graphviz(os, g_, boost::make_label_writer(boost::get(&VertexProperties::index,g_)), boost::make_label_writer(boost::get(&EdgeProperties::cap_info,g_)), boost::default_writer(), propmapIndex);
    //boost::write_graphviz(os, g_, boost::make_label_writer(boost::get(&VertexProperties::index,g_)), boost::make_label_writer(boost::get(&EdgeProperties::capacity,g_)), boost::default_writer(), propmapIndex);
    //boost::write_graphviz(os, g_, boost::make_label_writer(boost::get(&VertexProperties::index,g_)), boost::make_label_writer(boost::get(&EdgeProperties::residual_capacity,g_)), boost::default_writer(), propmapIndex);
    //boost::write_graphviz(os, g_, boost::make_label_writer(boost::get(&VertexProperties::index,g_)), boost::make_label_writer(boost::get(&EdgeProperties::v,g_)), boost::default_writer(), propmapIndex);
    //boost::write_graphviz(os, g_, boost::make_label_writer(boost::get(&VertexProperties::index,g_)), boost::make_label_writer(boost::get(&EdgeProperties::distance,g_)), boost::default_writer(), propmapIndex);
    //boost::write_graphviz(os, g_, boost::make_label_writer(boost::get(&VertexProperties::index,g_)), boost::make_label_writer(boost::get(&EdgeProperties::fatness,g_)), boost::default_writer(), propmapIndex);
}

void UnitigGraph::debug()
{
	// DEBUG
	auto&& numV = boost::num_vertices(g_);
	auto&& numE = boost::num_edges(g_);
	
	std::cerr << "Unitig graph has " << numV << " vertices and " << numE << " edges, starting cleaning" << std::endl;
}
