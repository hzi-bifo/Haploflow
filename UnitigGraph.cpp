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
UnitigGraph::UnitigGraph() : cc_(1), k_(0), read_length_(0)
{
}
// constructor of the so-called UnitigGraph
// unifies all simple paths in the deBruijnGraph to a single source->sink path
// all remaining nodes have either indegree != outdegree or indegree == outdegree > 1
UnitigGraph::UnitigGraph(deBruijnGraph& dbg, std::string p, float error_rate) : cc_(1), k_(dbg.k_), read_length_(dbg.read_length_)
{
	std::cerr << "deBruijnGraph has " << dbg.getSize() << " vertices" << std::endl;
	std::cerr << "Building unitig graph from deBruijn graph..." << std::endl;
    
	clock_t t = clock();
	auto&& junc = dbg.getJunctions();
	unsigned int index = 1;
	auto&& out_unbalanced = junc.first;
	auto&& in_unbalanced = junc.second;
    auto&& thresholds = calculate_thresholds(dbg, p, error_rate);
    thresholds_ = thresholds;
    graph_map_.resize(thresholds.size());
    graphs_.resize(thresholds.size());
    for (auto&& g : graphs_)
    {
        g = new UGraph;
    }
    // starting from the sources, we build the unitig graph
	for (auto& v : out_unbalanced)
	{
		std::string curr = v.get_kmer();
		Vertex* source = dbg.getVertex(curr);
        unsigned int cc = source->cc - 1; // cc start at 1
        float threshold = thresholds.at(cc);
		// make a guess whether we are relevant already 
		if (!source->is_flagged() and threshold > 0 and source->get_total_out_coverage() >= threshold)
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
        float threshold = thresholds.at(cc);
		if (!source->is_flagged() and threshold > 0 and source->get_total_in_coverage() >= threshold)
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
    unsigned int total_size = 0;
    for (unsigned int cc = 0; cc < graphs_.size(); cc++)
    {
        total_size += boost::num_vertices(*(graphs_.at(cc)));
    }
    std::cerr << "Unitig graph has " << total_size << " vertices" << std::endl;
}

UnitigGraph::~UnitigGraph()
{
    for (auto&& g : graphs_)
    {
        delete g;
    }
}

std::vector<float> UnitigGraph::calculate_thresholds(deBruijnGraph& dbg, std::string path, float error_rate)
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
    return get_thresholds(cov_distr, path, error_rate);
}

std::vector<float> UnitigGraph::get_thresholds(std::vector<std::map<unsigned int, unsigned int>>& cov_distr, std::string path, float error_rate)
{
    std::vector<float> thresholds;
    unsigned int i = 0;
    for (auto& covs : cov_distr)
    {
        if (covs.empty())
        {
            continue;
        }
        unsigned int members = 0;
        std::vector<float> sorted_coverage;
        sorted_coverage.resize(covs.rbegin()->first, 0.f); // get last (= biggest) element of map
        std::string filename = path + "Cov" + std::to_string(i) + ".tsv";
        std::ofstream outfile (filename);
        for (auto&& cov : covs)
        {
            auto pos = cov.first - 1;
            auto val = cov.second;
            members += val;
            sorted_coverage[pos] = val;
            outfile << pos << '\t' << val << std::endl;
        }
        if (members < 150) //less than 500 kmers
        {
            i++;
            thresholds.push_back(std::numeric_limits<float>::max()); // skip graph in creation
            continue;
        }
        auto diffs = finite_difference(sorted_coverage);
        float turning_point = diffs.first; 
        float inflexion_point = diffs.second; 
        if (inflexion_point > 1 and turning_point > 1)
        {
            std::cerr << "Graph " << i << " threshold set to " << std::max(1.f, std::min(turning_point, inflexion_point)) << std::endl;
        }
        /*if (std::min(inflexion_point, turning_point) > 1.f)
        {
            thresholds.push_back(2.f);
        }
        else
        {
            thresholds.push_back(1.f);
        }*/
        thresholds.push_back(std::max(1.f, std::min(turning_point, inflexion_point))); // the threshold should never be less than 1
        i++;
    }
    return thresholds;
}

// calculating finite differences, can only to 1st (false) and 2nd (true) order
std::pair<float, float> UnitigGraph::finite_difference(std::vector<float> input)
{
    if (input.size() < 3)
    {
        return std::make_pair(1.f, 1.f);
    }
    float inflexion1 = 0;
    float inflexion2 = 0;
    for (unsigned int i = 2; i < input.size(); i++)
    {
        float diff1 = 0.f;
        float diff2 = 0.f;
        diff1 = input[i - 1] - input[i]; //1st finite difference
        diff2 = input[i - 2] - 2 * input[i - 1] + input[i]; // 2nd finite derivative
        if (diff1 < 0 and inflexion1 == 0)
        {
            inflexion1 = float(i - 1);// first diff is shifted by 1
        }
        if (diff2 < 0 and inflexion2 == 0)
        {
            inflexion2 = float(i - 2);// second diff is shifted by 2
        }
    }
    if (inflexion1 == 0)
    {
        inflexion1 = 1.f;
    }
    if (inflexion2 == 0)
    {
        inflexion2 = 1.f;
    }
    return std::make_pair(inflexion1, inflexion2);
}

// calculates cumulative minimum of in
std::vector<float> UnitigGraph::cummin(std::vector<float>& in, unsigned int pos)
{
    std::vector<float> cummin;
    cummin.reserve(in.size());
    float curr = in[pos];
    for (unsigned int i = 0; i < pos; i++)
    {
        cummin.push_back(0);
    }
    for (auto it = in.begin() + pos; it != in.end(); ++it)
    {
        auto val = *it;
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
    UGraph* g = graphs_.at(ccc - 1);
	UVertex uv = boost::add_vertex(*g);
	(*index)++;
    // set vertex properties
    (*g)[uv].name = name;
	(*g)[uv].index = *index;
    (*g)[uv].scc = 1;
    (*g)[uv].tarjan_index = 0; // needs to be 0 to find out whether it has been set
    (*g)[uv].onStack = false;
    (*g)[uv].cc = ccc - 1; // set cc to current CC (startin*(g) with 0 here)
    (*g)[uv].visiting_time = 0;
	
    auto&& ins = std::make_pair(*index,uv);
	graph_map_.at(ccc - 1).insert(ins);
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
			uv = graph_map_.at(junction->cc - 1).at(idx);
		}
		junction->visit(); // make sure the next time we find it we dont add it another time
		auto&& following = addNeighbours(seq, succ, pred, dbg, index, uv, threshold, error_rate); // finding the next unbalanced vertices
		for (auto v : following)
		{	
			// if no sequence is returned, no vertex was added, so we do not need to continue on this vertex
			if (v.second != "")
				todo.push(v);
		}
	} while (!todo.empty());
}

// iterating over all neighbours of current node, build the different sequences to the next unbalanced node
std::vector<std::pair<Vertex*,std::string> > UnitigGraph::addNeighbours(std::string& curr, const std::vector<char>& succ, const std::vector<char>& pred, deBruijnGraph& dbg, unsigned int* index, UVertex& uv, float threshold, float error_rate)
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
                    curr_start = currV->get_read_starts() * currV->get_out_coverage(n)/total_out;
                next = curr.substr(1) + n;
                Vertex* nextV = dbg.getVertex(next);
                Sequence s = *dbg.getSequence(next);
                sequence += n;
                if (!nextV->is_flagged() and currV->get_out_coverage(n) > error_rate * total_out)
                {
                    following.push_back(buildEdge(uv, nextV, next, sequence, index, coverage, pcov, dbg, curr_start, threshold));
                }
            }
            // if we are a reverse complement, we actually want to add the path in reverse order
            else
            {
                float curr_end = 0.f;
                if (currV->get_total_in_coverage() != 0)
                    curr_end = currV->get_read_ends() * currV->get_in_coverage(deBruijnGraph::complement(n))/total_in;
                next = deBruijnGraph::complement(n) + curr.substr(0,curr.length() - 1);
                Vertex* nextV = dbg.getVertex(next);
                Sequence s = *dbg.getSequence(next);
                sequence += curr.back(); // the predecessor points to the current vertex with the last char of curr (by definition)
                if (!nextV->is_flagged() and currV->get_out_coverage(n) > error_rate * total_out)
                {
                    following.push_back(buildEdgeReverse(uv, nextV, next, sequence, index, coverage, pcov, dbg, curr_end, threshold));
                }
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
                    curr_end = currV->get_read_ends() * currV->get_in_coverage(n)/total_in;
                prev = n + curr.substr(0,curr.length() - 1);
                Vertex* nextV = dbg.getVertex(prev);
                Sequence s = *dbg.getSequence(prev);
                sequence += curr.back(); // the predecessor points to the current vertex with the last char of curr
                if (!nextV->is_flagged() and currV->get_in_coverage(n) > error_rate * total_in)
                {
                    following.push_back(buildEdgeReverse(uv, nextV, prev, sequence, index, coverage, pcov, dbg, curr_end, threshold));
                }
            }
            else
            {
                float curr_start = 0.f;
                if (currV->get_total_out_coverage() != 0)
                    curr_start = currV->get_read_starts() * currV->get_out_coverage(deBruijnGraph::complement(n))/total_out;
                prev = curr.substr(1) + deBruijnGraph::complement(n);
                Vertex* nextV = dbg.getVertex(prev);
                Sequence s = *dbg.getSequence(prev);
                sequence += deBruijnGraph::complement(n);
                if (!nextV->is_flagged() and currV->get_in_coverage(n) > error_rate * total_in)
                {
                    following.push_back(buildEdge(uv, nextV, prev, sequence, index, coverage, pcov, dbg, curr_start, threshold));
                }
            }
        }
    }
	currV->flag(); // this vertex is done
	return following;
}

// go back through the graph until the next unbalanced node is found and add an ("reversed") edge
std::pair<Vertex*,std::string> UnitigGraph::buildEdgeReverse(UVertex trg, Vertex* nextV, std::string prev, std::string& sequence, unsigned int* index, float coverage, float pcov, deBruijnGraph& dbg, float curr_end, float threshold)
{
    unsigned int cc = nextV->cc - 1;
    UGraph* g_ = graphs_.at(cc);
	char lastchar = (*g_)[trg].name.back(); //char with which we are pointing to trg
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
	if (!nextV->is_visited()) // TODO if coverage is low but the (unique) sequence is long, still add
	{// if the next vertex has been visited it already is part of the unitiggraph, otherwise add it
		nextV->visit();
		addVertex(index, prev, nextV->cc); // the vertex is new and found to be relevant
		nextV->index = *index;
	}
	else if (nextV->index == 0)
	{
		return std::make_pair(nextV,""); // path has too low coverage
	}
	UVertex src = graph_map_.at(nextV->cc - 1).at(nextV->index);
	auto e = boost::edge(src,trg,*g_);
	// if edge has been added or the immediate neighbour is an unbalanced vertex, do not add edge
	if ((!e.second or (e.second and (*g_)[e.first].name != sequence)))
	{
        //set new edge's information
		e = boost::add_edge(src,trg,*g_);
		std::reverse(sequence.begin(), sequence.end()); // we add the path from the found node to trg, the sequence was added in reverse order
		std::string old_name = (*g_)[e.first].name;
		if (e.second) // TODO v-S->w-T->v is treated like v<-S-w<-T-v (should be ST self-loop and TS self-loop)
		{
            (*g_)[e.first].name = old_name + sequence;
		}
		else
		{
            (*g_)[e.first].name = sequence;    
		}
        (*g_)[e.first].last_visit = 0;
        (*g_)[e.first].capacity = avg;
        (*g_)[e.first].residual_capacity = avg;
        (*g_)[e.first].cap_info.avg = avg;
        (*g_)[e.first].cap_info.max = max;
        (*g_)[e.first].cap_info.min = min;
        (*g_)[e.first].cap_info.first = first;
        (*g_)[e.first].cap_info.last = last;
        (*g_)[e.first].cap_info.length = (*g_)[e.first].name.length();
        (*g_)[e.first].residual_cap_info.avg = avg;
        (*g_)[e.first].residual_cap_info.max = max;
        (*g_)[e.first].residual_cap_info.min = min;
        (*g_)[e.first].residual_cap_info.first = first;
        (*g_)[e.first].residual_cap_info.last = last;
        (*g_)[e.first].residual_cap_info.length = (*g_)[e.first].name.length();
        (*g_)[e.first].visited = false;
	}
	return std::make_pair(nextV,prev);
}

// same function like the reverse one, but going forward and finding successors
std::pair<Vertex*,std::string> UnitigGraph::buildEdge(UVertex src, Vertex* nextV, std::string next, std::string& sequence, unsigned int* index, float coverage, float pcov, deBruijnGraph& dbg, float curr_start, float threshold)
{
	// with a little effort this can be moved inside the while loop for efficiency reasons
    unsigned int cc = nextV->cc - 1;
    UGraph* g_ = graphs_.at(cc);
    float first_char = (*g_)[src].name.front();
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
	if (!nextV->is_visited())
	{
		nextV->visit();
		addVertex(index, next, nextV->cc);
		nextV->index = *index;
	}
	else if (nextV->index == 0)
	{
		return std::make_pair(nextV,"");
	}
	UVertex trg = graph_map_.at(nextV->cc - 1).at(nextV->index);
	auto e = boost::edge(src,trg,*g_);
	if ((!e.second or (e.second and (*g_)[e.first].name != sequence)))
	{
		e = boost::add_edge(src,trg,*g_);
		auto old_name = (*g_)[e.first].name;
		if (e.second) // TODO v-S->w-T->v is treated like v<-S-w<-T-v (should be ST self-loop and TS self-loop)
		{
            (*g_)[e.first].name = sequence + old_name;
		}
		else
		{
            (*g_)[e.first].name = sequence;
		}
        (*g_)[e.first].last_visit = 0;
        (*g_)[e.first].capacity = avg;
        (*g_)[e.first].residual_capacity = avg;
        (*g_)[e.first].cap_info.avg = avg;
        (*g_)[e.first].cap_info.max = max;
        (*g_)[e.first].cap_info.min = min;
        (*g_)[e.first].cap_info.first = first;
        (*g_)[e.first].cap_info.last = last;
        (*g_)[e.first].cap_info.length = (*g_)[e.first].name.length();
        (*g_)[e.first].residual_cap_info.avg = avg;
        (*g_)[e.first].residual_cap_info.max = max;
        (*g_)[e.first].residual_cap_info.min = min;
        (*g_)[e.first].residual_cap_info.first = first;
        (*g_)[e.first].residual_cap_info.last = last;
        (*g_)[e.first].residual_cap_info.length = (*g_)[e.first].name.length();
        (*g_)[e.first].visited = false;
	}
	return std::make_pair(nextV,next);
}

// contracts all simple paths in graph to a single source-sink connection
void UnitigGraph::contractPaths(unsigned int cc)
{
    UGraph* g_ = graphs_.at(cc);
    boost::graph_traits<UGraph>::vertex_iterator vi, vi_end, next;
    boost::tie(vi, vi_end) = boost::vertices(*g_);
    for (next = vi; vi != vi_end; vi = next)
    {
        ++next;
        unsigned int indegree = boost::in_degree(*vi, *g_);
        unsigned int outdegree = boost::out_degree(*vi, *g_);

        // if in and outdegree is 1, we are on a simple path and can contract again
        if (outdegree == 1 and indegree == 1)
        {
            auto&& ie = boost::in_edges(*vi,*g_);
            auto&& oe = boost::out_edges(*vi,*g_);
            auto&& new_source = boost::source(*ie.first,*g_);
            auto&& new_target = boost::target(*oe.first,*g_);
            if (new_source == new_target)
            {
                continue; // do not contract to single vertex (which might get deleted)
            }
            auto&& e = boost::edge(new_source,*vi,*g_);
            auto&& f = boost::edge(*vi, new_target,*g_); // coverage etc of second edge to be contracted
            std::string seq = (*g_)[e.first].name;
            unsigned int w = seq.length();

            Capacity cap_info_e = (*g_)[e.first].cap_info;
            Capacity cap_info_f = (*g_)[f.first].cap_info;
            float max = std::max(cap_info_e.max, cap_info_f.max);
            float min = std::min(cap_info_e.min, cap_info_f.min);
            float first = cap_info_e.first;
            float last = cap_info_f.last;

            // TODO temporarily disabled to study effect of not contracting
            //if (std::abs(cap_info_f.first - cap_info_e.last) > threshold_ or std::abs(cap_info_f.first - cap_info_e.last) > threshold_) // do not contract paths which have high divergence in capacity
            //    continue;
            
            float capacity = (*g_)[e.first].capacity * w;
            e = boost::edge(*vi,new_target,*g_);
            seq += (*g_)[e.first].name; // append the sequence
            capacity += (*g_)[e.first].capacity * (seq.length() - w);
            if (seq.length() > 0)
                capacity /= seq.length(); // currently using the average coverage on the contracted path

            auto&& new_e = boost::add_edge(new_source,new_target,*g_);
            (*g_)[new_e.first].last_visit = 0;
            (*g_)[new_e.first].name = seq;
            (*g_)[new_e.first].capacity = capacity;
            (*g_)[new_e.first].residual_capacity = capacity;
            (*g_)[new_e.first].cap_info.avg = capacity;
            (*g_)[new_e.first].cap_info.max = max;
            (*g_)[new_e.first].cap_info.min = min;
            (*g_)[new_e.first].cap_info.first = first;
            (*g_)[new_e.first].cap_info.last = last;
            (*g_)[new_e.first].cap_info.length = (*g_)[new_e.first].name.length();
            (*g_)[new_e.first].residual_cap_info.avg = capacity;
            (*g_)[new_e.first].residual_cap_info.max = max;
            (*g_)[new_e.first].residual_cap_info.min = min;
            (*g_)[new_e.first].residual_cap_info.first = first;
            (*g_)[new_e.first].residual_cap_info.last = last;
            (*g_)[new_e.first].residual_cap_info.length = (*g_)[e.first].name.length();
            (*g_)[new_e.first].visited = false;
            boost::clear_vertex(*vi,*g_);
            boost::remove_vertex(*vi,*g_);
        }
    }
}

void UnitigGraph::removeStableSets(unsigned int cc)
{
    UGraph* g_ = graphs_.at(cc);
	boost::graph_traits<UGraph>::vertex_iterator vi, vi_end, next;
	boost::tie(vi, vi_end) = boost::vertices(*g_);
	for (next = vi; vi != vi_end; vi = next) {
		++next;
		unsigned int indegree = boost::in_degree(*vi, *g_);
		unsigned int outdegree = boost::out_degree(*vi,*g_);
		if (outdegree == 0 and indegree == 0)
		{
			boost::remove_vertex(*vi,*g_);
		}
	}
}


bool UnitigGraph::hasRelevance(unsigned int cc)
{
    UGraph* g_ = graphs_.at(cc);
    unsigned int length = 0;
    for (auto&& e : boost::edges(*g_))
    {
        length += (*g_)[e].name.size();
    }
    if (length <= 150)
    {
        std::cerr << "Graph undercutting threshold of 150 characters (" << length << ")" << std::endl;
    }
    return length > 150;
}

// the graph might contain some unconnected vertices, clean up
void UnitigGraph::cleanGraph(unsigned int cc, float error_rate)
{
    UGraph* g_ = graphs_.at(cc);
    for (auto e : boost::edges(*g_))
    {
        (*g_)[e].last_visit = 0; // reusing for number of allowed paths
    }
    removeLow_cutEnds(cc, error_rate);
	removeStableSets(cc);
    contractPaths(cc);
    removeShortPaths(cc);
	removeStableSets(cc);
    contractPaths(cc);
    //expandCycles(cc);
}

void UnitigGraph::removeLow_cutEnds(unsigned int cc, float error_rate)
{
    UGraph* g_ = graphs_.at(cc);
    std::set<UEdge> toDelete;
    for (auto&& v : boost::vertices(*g_))
    {
        float out_degree = 0.f;
        float in_degree = 0.f;
        for (auto&& oe : boost::out_edges(v, *g_))
        {
            out_degree += (*g_)[oe].capacity;
        }
        for (auto&& ie : boost::in_edges(v, *g_))
        {
            in_degree += (*g_)[ie].capacity;
        }
        for (auto&& oe : boost::out_edges(v, *g_))
        {
            if ((*g_)[oe].cap_info.max < thresholds_[cc] or (*g_)[oe].capacity < error_rate * out_degree)
            {
                toDelete.insert(oe);
            }
        }
        for (auto&& ie : boost::in_edges(v, *g_))
        {
            if ((*g_)[ie].cap_info.max < thresholds_[cc] or (*g_)[ie].capacity < error_rate * in_degree)
            {
                toDelete.insert(ie);
            }
        }
    }
    for (auto&& e : toDelete)
    {
        boost::remove_edge(e, *g_);
    }
    toDelete.clear();
    contractPaths(cc);
    for (auto&& v : boost::vertices(*g_))
    {
        if (boost::in_degree(v, *g_) == 0)
        {
            for (auto&& oe : boost::out_edges(v, *g_))
            {
                if ((*g_)[oe].capacity < thresholds_[cc])
                {
                    toDelete.insert(oe);
                }
            }
        }
        else if (boost::out_degree(v, *g_) == 0)
        {
            for (auto&& ie : boost::in_edges(v, *g_))
            {
                if ((*g_)[ie].capacity < thresholds_[cc])
                {
                    toDelete.insert(ie);
                }
            }
        }
    }
}

/*void UnitigGraph::expandCycles(unsigned int cc)
{
    
}*/

void UnitigGraph::removeShortPaths(unsigned int cc)
{
    UGraph* g_ = graphs_.at(cc);
    std::vector<UEdge> toDelete;
    for (auto&& e : boost::edges(*g_))
    {
        auto source = boost::source(e, *g_);
        auto target = boost::target(e, *g_);
        // we are a simple edge between two vertices, unconnected to the rest
        if ((boost::out_degree(target, *g_) == 0 and boost::in_degree(source, *g_) == 0
            and boost::out_degree(source, *g_) == 1 and boost::in_degree(target, *g_) == 1))
        {
            if ((*g_)[e].name.size() < read_length_)
            {
                toDelete.push_back(e);
            }
        }
        else if ((boost::out_degree(target, *g_) == 1 and boost::in_degree(source, *g_) == 1
            and boost::out_degree(source, *g_) == 1 and boost::in_degree(target, *g_) == 1))
        {
            float length = (*g_)[e].name.size();
            auto rev_e = boost::edge(target, source, *g_);
            if (rev_e.second)
            {
                length += (*g_)[rev_e.first].name.size();
                if (length < read_length_)
                {
                    toDelete.push_back(e);
                }
            }
        }
    }
    for (auto&& e : toDelete)
    {
        boost::remove_edge(e, *g_);
    }
}

// Tests whether two percentages "belong together"
bool UnitigGraph::test_hypothesis(float to_test_num, float to_test_denom, float h0, float threshold)
{
    float diff = std::abs(to_test_num - h0 * to_test_denom); // the absolute difference between num and expected num
    //bool less = to_test_num < to_test_denom * h0;
    //float perc = to_test_num/to_test_denom; //currently unused
    return (/*less or*/ diff < threshold); // they differ by less than threshold 
}


float UnitigGraph::in_capacity(UVertex source, unsigned int cc)
{
    UGraph* g_ = graphs_.at(cc);
    float capacity = 0;
    for (auto ie : boost::in_edges(source, *g_))
    {
        capacity += (*g_)[ie].capacity;
    }
    return capacity;
}

float UnitigGraph::out_capacity(UVertex target, unsigned int cc)
{
    UGraph* g_ = graphs_.at(cc);
    float capacity = 0;
    for (auto oe : boost::out_edges(target, *g_))
    {
        capacity += (*g_)[oe].capacity;
    }
    return capacity;
}

// run dijsktra with fatness as optimality criterion, marks the graph with the distances from seed
void UnitigGraph::dijkstra(UEdge seed, bool init, bool local, unsigned int cc)
{
    UGraph* g_ = graphs_.at(cc);
    auto edge_compare = [&](UEdge e1, UEdge e2){ //sort by biggest fatness
        return (*g_)[e1].fatness < (*g_)[e2].fatness;
    };
    std::vector<UEdge> q;
    if (!local)
    {
        for (auto e : boost::edges(*g_)) //initialise distances and fatness
        {
            (*g_)[e].fatness = 0;
            (*g_)[e].distance = std::numeric_limits<unsigned int>::max();
            if (e == seed) //
            {
                (*g_)[e].distance = 0;
                (*g_)[e].fatness = std::numeric_limits<float>::max();
                // such that source plays a lesser role
            }
        }
    }
    if (local) // if fatness is not initialised for all vertices, initialize for all following
    {
        q.push_back(seed);
        while (!q.empty()) // classic dijsktra routine (cancelling when cycle found)
        {
            auto curr = q.back();
            q.pop_back();
            (*g_)[curr].visited = true;
            (*g_)[curr].fatness = 0;
            auto target = boost::target(curr, *g_);
            for (auto oe : boost::out_edges(target, *g_))
            {
                if (!(*g_)[oe].visited)
                {
                    q.push_back(oe);
                }
            }
        }
    }
    (*g_)[seed].fatness = std::numeric_limits<float>::max();
    if ((*g_)[seed].distance == std::numeric_limits<float>::max())
    {
        (*g_)[seed].distance = 0;
    }
    for (auto e : boost::edges(*g_)) //unvisit for next run
    {
        (*g_)[e].visited = false;
    }
    q.push_back(seed);
    while (!q.empty()) // classic dijsktra routine (cancelling when cycle found)
    {
        auto curr = q.back();
        q.pop_back();
        //(*g_)[curr].visited = true;
        auto target = boost::target(curr, *g_);
        for (auto oe : boost::out_edges(target, *g_))
        {
            float fat = (*g_)[oe].fatness;
            if (fat < std::min((*g_)[curr].fatness, (init ? (*g_)[oe].residual_capacity : (*g_)[oe].capacity)))
            {
                (*g_)[oe].fatness = std::min((*g_)[curr].fatness, (init ? (*g_)[oe].residual_capacity : (*g_)[oe].capacity));
                (*g_)[oe].prev = curr;
                (*g_)[oe].distance = (*g_)[curr].distance + (*g_)[oe].name.size();
                auto old_pos = std::find(q.begin(), q.end(), oe);
                if (old_pos != q.end())
                {
                    q.erase(old_pos);
                }
                auto pos = std::lower_bound(q.begin(), q.end(), oe, edge_compare);
                q.insert(pos, oe);
            }
        }
    }
    (*g_)[seed].fatness = (init ? (*g_)[seed].residual_capacity : (*g_)[seed].capacity);
}

std::pair<UEdge, float> UnitigGraph::get_target(UEdge seed, bool lenient, unsigned int cc)
{
    UGraph* g_ = graphs_.at(cc);
    float max_dist = 0;
    float running_distance = 0;
    UEdge last;
    auto visits = (*g_)[seed].visits;
    unvisit(cc);
    std::vector<UEdge> q = {seed};
    while (!q.empty()) // classic dijsktra routine (cancelling when cycle found)
    {
        auto curr = q.back();
        q.pop_back();
        
        bool same_visit = false;
        
        for (auto v : (*g_)[curr].visits)
        {
            same_visit = (std::find(visits.begin(), visits.end(), v) != visits.end());
            if (same_visit)
                break;
        }

        (*g_)[curr].visited = true;
        auto target = boost::target(curr, *g_);
        /*if (same_visit and num_visits == 1 and (*g_)[curr].distance < std::numeric_limits<unsigned int>::max() and (running_distance == 0 or ((*g_)[curr].fatness > 0 and running_fatness/(*g_)[curr].fatness < (*g_)[curr].distance/running_distance)))
        {
            running_fatness = (*g_)[curr].fatness;
            running_distance = (*g_)[curr].distance;
            last = curr;
            max_dist = (*g_)[curr].distance;
        }*/
        if ((*g_)[curr].distance > running_distance and (*g_)[curr].distance != std::numeric_limits<unsigned int>::max() and same_visit)
        {
            running_distance = (*g_)[curr].distance;
            last = curr;
            max_dist = (*g_)[curr].distance;
        }
        for (auto oe : boost::out_edges(target, *g_))
        {
            if (!(*g_)[oe].visited and (*g_)[oe].prev == curr)
            {
                q.push_back(oe);
            }
        }
    }
    if (max_dist == 0)
    {
        last = seed; // so we don't return nothing
    }
    unvisit(cc);
    return std::make_pair(last, max_dist);
}

std::vector<UEdge> UnitigGraph::fixFlow(UEdge seed, unsigned int cc)
{
    unvisit(cc);
    dijkstra(seed, false, false, cc);
    auto path = find_fattest_path(seed, cc);
    UGraph* g_ = graphs_.at(cc);
    bool corrected = true;
    unsigned int pos = 0;
    bool eq = false;
    std::vector<UEdge> return_path;
    while (corrected and !eq)
    {
        std::vector<UEdge> tmp_path;
        bool dip = false;
        unsigned int i = 0;
        for (auto& e : path)
        {
            tmp_path.push_back(e);
            if (i > pos and i < path.size() - 1)
            {
                auto trg = boost::target(e, *g_);
                bool less = false;
                bool more = false;
                UEdge next = path[i + 1];
                for (auto&& oe : boost::out_edges(trg, *g_))
                {
                    if (oe != next)
                    {
                        dip = dip or ((*g_)[oe].fatness == (*g_)[next].fatness);
                        less = (*g_)[next].capacity < 1.05 * (*g_)[oe].capacity or (*g_)[next].capacity < (*g_)[oe].capacity + thresholds_[cc];
                        more = (*g_)[next].capacity * 1.05 > (*g_)[oe].capacity or (*g_)[next].capacity + thresholds_[cc] > (*g_)[oe].capacity;
                        eq = less and more;
                    }
                    if (eq)
                    {
                        break;
                    }
                }
                if (eq) // if two out edges are the same capacity: break immediately
                {
                    break;
                }
                if (dip) // if dip in between caused problems: try to fix flow downstream
                {
                    unvisit(cc);
                    pos = i;
                    dijkstra(e, false, true, cc);
                    break;
                }
            }
            i++;
            return_path.push_back(e);
        }
        return_path = tmp_path;
        corrected = (dip and i == pos);
    }
    if (!eq)
    {
        unvisit(cc);
        return_path = find_fattest_path(seed, cc);
    }
    return return_path;
}

// Calculates the fattest path through the graph, marks vertices as being visited
std::vector<UEdge> UnitigGraph::find_fattest_path(UEdge seed, unsigned int cc)
{
    UGraph* g_ = graphs_.at(cc);
    //auto source = boost::source(seed, *g_);
    //auto target = boost::target(seed, *g_);
    //std::cerr << "Source: " << (*g_)[source].index << " -> " << (*g_)[target].index << ": " << (*g_)[seed].capacity << std::endl;
    auto trg = get_target(seed, false, cc);
    auto last = trg.first;
    //auto source2 = boost::source(last, *g_);
    //auto target2 = boost::target(last, *g_);
    //std::cerr << "Target: " << (*g_)[source2].index << " -> " << (*g_)[target2].index << ": " << (*g_)[last].capacity << std::endl;
    float max_dist = trg.second;
    if (max_dist == 0) // path is only one edge, no longest path
    {
        return std::vector<UEdge>{seed};
    }
    auto curr = last;
    std::deque<UEdge> path = {curr};
    (*g_)[curr].last_visit++;
    while (curr != seed and (*g_)[curr].distance > 0 and (*g_)[curr].distance < std::numeric_limits<unsigned int>::max() and !(*g_)[curr].visited) // this means the distance has been set, i.e. the vertex has been reached
    {
        (*g_)[curr].visited = true;
        curr = (*g_)[curr].prev;
        (*g_)[curr].last_visit++;
        path.push_front(curr);
    }
    //source = boost::source(path.front(), *g_);
    unsigned int i = 0;
    unsigned int j = 1;
    float seq_length = 0;
    std::vector<UEdge> ret;
    
    //std::cerr << (*g_)[source].index; 
    for (auto e : path)
    {
        //auto trg = boost::target(e, *g_);
        //std::cerr << " -> " << (*g_)[trg].index << "(" << (*g_)[e].capacity << " " << (*g_)[e].fatness << ") ";
        seq_length += (*g_)[e].name.size();
        auto ct = 0;
        auto avg = 0.f;
        for (unsigned int k = i; k < j; k++)
        {
            unsigned int visits = (*g_)[path[k]].visits.size();
            if (visits < 2)
            {
                avg += (*g_)[path[k]].capacity;
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
            seq_length -= (*g_)[path[i]].name.size();
            i++;
        }
        j++;
        ret.push_back(e);
    }
    //std::cerr << std::endl;
    return ret;
}

void UnitigGraph::unvisit(unsigned int cc)
{
    UGraph* g_ = graphs_.at(cc);
    for (auto e : boost::edges(*g_))
    {
        (*g_)[e].visited = false;
    }
}

float UnitigGraph::reduce_flow(std::vector<UEdge>& path, std::set<unsigned int>& unique_paths, unsigned int cc, bool init)
{
    UGraph* g_ = graphs_.at(cc);
    float removed_coverage = (*g_)[path.front()].capacity;
    float average = 0.f;
    unsigned int len = 0;
    for (auto e : path)
    {
        if (!init)
        {
            // first find out which path we are on (we delete this because it has been used then)
            auto to_remove = (*g_)[e].visits.begin();
            for (auto p : unique_paths)
            {
                to_remove = std::find((*g_)[e].visits.begin(), (*g_)[e].visits.end(), p);
                if (to_remove != (*g_)[e].visits.end())
                {
                    break; // TODO
                }
            }
            // TODO fix multiple "unique" paths
            if (to_remove != (*g_)[e].visits.end())
            {
                (*g_)[e].visits.erase(to_remove);
            }
            else
            {
                if (!(*g_)[e].visits.empty())
                {
                    (*g_)[e].visits.erase((*g_)[e].visits.begin());
                }
            }
        }
        // now check whether this was the last path or there are remaining paths to reduce capacity
        float val = init ? (*g_)[e].residual_capacity : (*g_)[e].capacity;
        float val_last = init ? (*g_)[e].residual_cap_info.last : (*g_)[e].cap_info.last;
        float val_first = init ? (*g_)[e].residual_cap_info.first : (*g_)[e].cap_info.first;

        auto threshold = thresholds_.at(cc);
        
        bool not_decreasing = val_first <= 1.1 * val_last or std::abs(val_first - val_last) < threshold; //test_hypothesis((*g_)[e].cap_info.first, (*g_)[e].cap_info.last, 1.2, threshold_); //not decreasing (first/last >= 1.2)
        bool not_increasing = val_last <= 1.1 * val_first or std::abs(val_first - val_last) < threshold; //test_hypothesis((*g_)[e].cap_info.last, (*g_)[e].cap_info.first, 1.2, threshold_); // not increasing (last/first >= 1.2)
        //cannot both be false, but both be true (if close to 1.2)
        
        if (val > threshold and not_decreasing and not_increasing)
        {
            if (!init)
            {
                if ((*g_)[e].visits.empty() or (*g_)[e].capacity <= threshold)
                {
                    (*g_)[e].capacity = 0;
                    (*g_)[e].cap_info.first = 0;
                    (*g_)[e].cap_info.last = 0;
                }
                else
                {
                    (*g_)[e].capacity = std::max(threshold, (*g_)[e].capacity - removed_coverage); // there might be paths, so leave a small amount
                    (*g_)[e].cap_info.first = std::max(threshold, (*g_)[e].cap_info.first - removed_coverage);
                    (*g_)[e].cap_info.last = std::max(threshold, (*g_)[e].cap_info.last - removed_coverage);
                }
                (*g_)[e].cap_info.avg = (*g_)[e].capacity;
            }
            else
            {
                (*g_)[e].residual_capacity = std::max(threshold, (*g_)[e].residual_capacity - removed_coverage); // there might be paths, so leave a small amount
                (*g_)[e].residual_cap_info.first = std::max(threshold, (*g_)[e].residual_cap_info.first - removed_coverage);
                (*g_)[e].residual_cap_info.last = std::max(threshold, (*g_)[e].residual_cap_info.last - removed_coverage);
                (*g_)[e].residual_cap_info.avg = (*g_)[e].residual_capacity;
            }
            removed_coverage = val - (init ? (*g_)[e].residual_capacity : (*g_)[e].capacity);
        }
        else if (val > threshold and not_decreasing and !not_increasing)
        {
            if (!init)
            {
                if ((*g_)[e].visits.empty() or (*g_)[e].capacity <= threshold)
                {
                    (*g_)[e].capacity = 0;
                    (*g_)[e].cap_info.first = 0;
                    (*g_)[e].cap_info.last = 0;
                }
                else
                {
                    (*g_)[e].capacity = std::max(threshold, (*g_)[e].cap_info.first - removed_coverage);
                    (*g_)[e].cap_info.last = std::max(threshold, (*g_)[e].cap_info.last - (val - (*g_)[e].capacity));
                }
                (*g_)[e].cap_info.avg = (*g_)[e].capacity;
                (*g_)[e].cap_info.first = (*g_)[e].capacity;
            }
            else
            {
                (*g_)[e].residual_capacity = std::max(threshold, (*g_)[e].residual_cap_info.first - removed_coverage);
                (*g_)[e].residual_cap_info.last = std::max(threshold, (*g_)[e].residual_cap_info.last - (val - (*g_)[e].residual_capacity));
                (*g_)[e].residual_cap_info.avg = (*g_)[e].residual_capacity;
                (*g_)[e].residual_cap_info.first = (*g_)[e].residual_capacity;
            }
            removed_coverage = val_last - (init ? (*g_)[e].residual_capacity : (*g_)[e].capacity);
        }
        else if (val > threshold and !not_decreasing and not_increasing)
        {
            if (!init)
            {
                if ((*g_)[e].visits.empty() or (*g_)[e].capacity <= threshold)
                {
                    (*g_)[e].capacity = 0;
                    (*g_)[e].cap_info.first = 0;
                    (*g_)[e].cap_info.last = 0;
                }
                else
                {
                    (*g_)[e].capacity = std::max(threshold, (*g_)[e].cap_info.last - removed_coverage);
                    (*g_)[e].cap_info.first = std::max(threshold, (*g_)[e].cap_info.first - (val - (*g_)[e].capacity));
                }
                (*g_)[e].cap_info.last = (*g_)[e].capacity;
                (*g_)[e].cap_info.avg = (*g_)[e].capacity;
            }
            else
            {
                (*g_)[e].residual_capacity = std::max(threshold, (*g_)[e].residual_cap_info.last - removed_coverage);
                (*g_)[e].residual_cap_info.first = std::max(threshold, (*g_)[e].residual_cap_info.first - (val - (*g_)[e].residual_capacity));
                (*g_)[e].residual_cap_info.last = (*g_)[e].residual_capacity;
                (*g_)[e].residual_cap_info.avg = (*g_)[e].residual_capacity;
            }
            removed_coverage = val_first - (init ? (*g_)[e].residual_capacity : (*g_)[e].capacity);
        }
        else
        {
            if (!init)
            {
                (*g_)[e].cap_info.first = 0;
                (*g_)[e].cap_info.last = 0;
                (*g_)[e].capacity = 0;
                (*g_)[e].cap_info.avg = 0;
            }
            else
            {
                (*g_)[e].residual_cap_info.first = 0;
                (*g_)[e].residual_cap_info.last = 0;
                (*g_)[e].residual_capacity = 0;
                (*g_)[e].residual_cap_info.avg = 0;
            }
            removed_coverage = val;
        }
        average *= len;
        average += removed_coverage * (*g_)[e].name.size();
        len += (*g_)[e].name.size();
        average /= len;
        if (init)
        {
            removed_coverage = std::max(removed_coverage, average);
        }
    }
    return average;
}

// Given all the chosen edges and their coverage fraction, builds the contigs and reduces flow accordingly
std::pair<std::string, float> UnitigGraph::calculate_contigs(std::vector<UEdge>& path, std::vector<float>& flows, unsigned int cc)
{
    UGraph* g_ = graphs_.at(cc);
    float flow = 0.;
    float max_flow = 0.;
    unsigned int i = 0;
    UEdge curr = path.front();
    UVertex source = boost::source(curr, *g_);
    std::set<unsigned int> paths;
    std::string contig = (*g_)[source].name;
    std::cerr << (*g_)[source].index;
    for (auto e : path)
    {
        contig += (*g_)[e].name;
        UVertex target = boost::target(e, *g_);
        std::cerr << " -> " << (*g_)[target].index << " (" << (*g_)[e].capacity << ", " << contig.size() << ")";
        if ((*g_)[e].visits.size() == 1)
        {
            if ((*g_)[e].capacity > max_flow)
            {
                max_flow = (*g_)[e].capacity;
            }
            paths.insert((*g_)[e].visits[0]);
            flow += (*g_)[e].capacity;
            i++;
        }
    }
    std::cerr << std::endl;
    if (i > 0)
    {
        flow /= i; // average flow over unqiue edges of path
    }
    else // no unique edges on path, remove duplicate path as far as possible and calculate new flow based on removed path
    {
        flow = 0;
        std::vector<unsigned int> to_remove = (*g_)[path.front()].visits; // choose one of these visits to be removed;
        for (auto e : path)
        {
            if ((*g_)[e].visits.size() > 1)
            {
                bool removed = false;
                for (auto v : to_remove)
                {
                    auto visit = std::find((*g_)[e].visits.begin(), (*g_)[e].visits.end(), v);
                    if (visit != (*g_)[e].visits.end())
                    {
                        (*g_)[e].visits.erase(visit);
                        removed = true;
                        flow += flows[*visit - 1];
                        break;
                    }
                }
                if (!removed)
                {
                    flow += flows[(*g_)[e].visits.front() - 1];
                    (*g_)[e].visits.erase((*g_)[e].visits.begin());
                }
            }
        }
        std::cerr << "Removed duplicated contig " << std::endl;
        std::cerr << contig << std::endl;
        contig = "";
        if (path.size() > 0)
            flow /= path.size();
    }
    reduce_flow(path, paths, cc, false);
    return std::make_pair(contig, flow);
}

std::vector<UEdge> UnitigGraph::blockPath(UEdge curr, unsigned int visits, unsigned int cc)
{
    UGraph* g_ = graphs_.at(cc);
    std::vector<UEdge> blockedPath;
    while (true)
    {
        (*g_)[curr].visits.push_back(visits);
        (*g_)[curr].last_visit = visits;
        blockedPath.push_back(curr);
        auto target = boost::target(curr, *g_);
        auto out_edges = boost::out_edges(target, *g_);
        float max = -1;
        float max_unvisited = -1;
        UEdge max_e;
        UEdge max_e_unvisited;
        for (auto&& e : out_edges)
        {
            if ((*g_)[e].visits.size() == 0)
            {
                if ((*g_)[e].residual_capacity > max_unvisited)
                {
                    max_unvisited = (*g_)[e].residual_capacity;
                    max_e_unvisited = e;
                }
            }
            if ((*g_)[e].residual_capacity > max) // maximal visited edge
            {
                max = (*g_)[e].residual_capacity;
                max_e = e;
            }
        }
        if (max != -1 and (*g_)[max_e].last_visit == visits) // the max edge has been visited in the same run
        { // TODO we still might want to continue, if the next edge has a similar coverage as the first had
            return blockedPath;
        }
        if (max_unvisited != -1)
        {
            curr = max_e_unvisited;
        }
        if (max != -1) // all out_edges are visited
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

UEdge UnitigGraph::get_next_source(unsigned int cc) /// just returns the highest capacity edge (TODO?)
{
    UGraph* g_ = graphs_.at(cc);
    auto edge_compare = [&](UEdge e1, UEdge e2){ //sort by biggest capacity
        return (*g_)[e1].capacity > (*g_)[e2].capacity;
    };
    UEdge source;

    auto sources = get_sources(cc);
    if (sources.size() > 0) // if there are sources, take highest possible source
    {
        std::sort(sources.begin(), sources.end(), edge_compare);
        source = sources.front();
    }
    else // else take highest unvisited edge
    {
        float max = 0;
        for (auto e : boost::edges(*g_))
        {
            if ((*g_)[e].capacity > max)
            {
                max = (*g_)[e].capacity;
                source = e;
            }
        }
    }
    return source;
}

std::pair<UEdge, bool> UnitigGraph::checkUnvisitedEdges(UEdge current_source, unsigned int cc)
{
    UGraph* g_ = graphs_.at(cc);
    bool unblocked = false;
    UEdge curr = current_source;
    auto source = boost::source(current_source, *g_);
    auto out_edges = boost::out_edges(source, *g_);
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
        if ((*g_)[next].visited) // to prevent cycling/multiple searches
            continue;
        else
            (*g_)[next].visited = true;
        if ((*g_)[next].last_visit == 0) // next edge has not been visited
        {
            if ((*g_)[next].residual_capacity > capacity) // and has higher capacity than other close-by edges
            {
                unblocked = true;
                curr = next;
                capacity = (*g_)[next].residual_capacity; // choose as next
            }
        }
        else // edge is visisted, continue search
        {
            source = boost::target(next, *g_);
            out_edges = boost::out_edges(source, *g_);
            for (auto e : out_edges)
            {
                if (!unblocked or (*g_)[e].last_visit == 0) // if we havent found an edge or next edge is unvisited: continue search
                    to_check.push(e);
            }
        }
    }
    return std::make_pair(curr, unblocked);
}

std::pair<UEdge, bool> UnitigGraph::getUnvisitedEdge(const std::vector<UEdge>& sources, unsigned int cc)
{
    UGraph* g_ = graphs_.at(cc);
    UEdge curr;
    bool unblocked = false;
    if (sources.size() == 0) // there is no source -> we are complete cycle -> pick highest capacity edge
    {
        float max = -1;
        for (auto e : boost::edges(*g_))
        {
            if ((*g_)[e].last_visit == 0 and (*g_)[e].residual_capacity > max)
            {
                unblocked = true;
                curr = e;
                max = (*g_)[e].residual_capacity;
            }
        }
    }
    else // at least one source which has been visited before
    {
        curr = sources.back();
        unblocked = ((*g_)[curr].last_visit == 0);
        for (auto e : sources) // check all sources for unchecked edges
        {
            auto nextUnvisited = checkUnvisitedEdges(e, cc);
            auto potential_source = nextUnvisited.first;
            if (nextUnvisited.second and ((*g_)[potential_source].residual_capacity > (*g_)[curr].residual_capacity or !unblocked))
            { // and check the highest capacity one
                unblocked = true;
                curr = potential_source;
            }
        }
    }
    return std::make_pair(curr, unblocked);
}

std::vector<UEdge> UnitigGraph::get_sources(unsigned int cc)
{
    UGraph* g_ = graphs_.at(cc);
    std::set<UEdge> sources;
    for (auto e : boost::edges(*g_))
    {
        auto src = boost::source(e, *g_);
        auto target = boost::target(e, *g_);
        auto in_degree = boost::in_degree(src, *g_);
        if (in_degree == 0) // only add all source edges if they havent beed added before
        {
            for (auto f : boost::out_edges(src, *g_))
            {
                sources.insert(f); // for every outedge of source start search
            }
        }
        else if (in_degree == 1 and src == target) // also add sources which have a self-loop
        {
            sources.insert(e);
        }
        bool eq1 = false;
        bool eq2 = false;
        if (in_degree == 1)
        {
            for (auto&& ie : boost::in_edges(src, *g_))
            {
                if (boost::source(ie, *g_) == target)
                {
                    eq1 = true;
                }
            }
            for (auto&& oe : boost::out_edges(target, *g_))
            {
                if (boost::target(oe, *g_) == src and boost::in_degree(target, *g_) == 1)
                {
                    eq2 = true;
                }
            }
            if (eq1 and eq2)
            {
                sources.insert(e);
            }
        }
    }
    return std::vector<UEdge>(sources.begin(), sources.end());
}

std::pair<std::vector<UEdge>, std::vector<float>> UnitigGraph::find_paths(unsigned int cc)
{
    UGraph* g_ = graphs_.at(cc);
    std::vector<UEdge> sources = get_sources(cc); //get sources of the graph (indegree = 0)
    auto edge_compare = [&](UEdge e1, UEdge e2){ //sort by biggest capacity
        return (*g_)[e1].capacity > (*g_)[e2].capacity;
    };
    std::sort(sources.begin(), sources.end(), edge_compare); // so that we search the highest source first
    unsigned int visits = 1;
    std::vector<std::vector<UEdge>> unique;
    std::vector<UEdge> started_from;
    std::vector<float> unique_paths;
    while (true)
    {
        unvisit(cc); // TODO only current cc
        std::pair<UEdge, bool> unvisited = getUnvisitedEdge(sources, cc);
        UEdge curr = unvisited.first;
        bool unblocked = unvisited.second;
        std::vector<UEdge> blockedPath;
        if (unblocked)
        {
            unique.push_back(std::vector<UEdge>{});
            dijkstra(curr, true, false, cc);
            started_from.push_back(curr); // add the edge from which we started
            blockedPath = blockPath(curr, visits, cc); //marks the first path
        }
        else
        {
            break;
        }
        unsigned int length = 0;
        for (auto e : boost::edges(*g_))
        {
            if (!(*g_)[e].visits.empty() and (*g_)[e].visits.front() == visits)
            {
                unique[visits - 1].push_back(e);
            }
        }
        // we now have the tentative paths, now check how many edges are unique per path
        unique_paths.push_back(remove_non_unique_paths(unique, blockedPath, length, visits - 1, cc));
        visits++;
    }
    for (auto&& e : boost::edges(*g_))
    {
        if ((*g_)[e].visits.empty())
        {
            (*g_)[e].visits.push_back(visits);
            started_from.push_back(e);
            unique_paths.push_back(0.);
            (*g_)[e].residual_capacity = 0; // there might be paths, so leave a small amount
            (*g_)[e].residual_cap_info.first = 0;
            (*g_)[e].residual_cap_info.last = 0;
            (*g_)[e].residual_cap_info.avg = 0;
        }
    }
    return std::make_pair(started_from, unique_paths);
}

float UnitigGraph::remove_non_unique_paths(std::vector<std::vector<UEdge>>& unique, std::vector<UEdge>& blockedPath, unsigned int length, unsigned int visits, unsigned int cc)
{
    UGraph* g_ = graphs_.at(cc);
    auto size = unique[visits].size();
    float median = 0.f;
    if (size < 0.02 * boost::num_edges(*g_) and size < 15 and length < 500) //TODO parameters
    {
        for (auto e : blockedPath)
        {
            if ((*g_)[e].visits.size() > 1) // so we dont not visit some edges
            {
                auto pos = std::find((*g_)[e].visits.begin(), (*g_)[e].visits.end(), visits + 1);
                if (pos != (*g_)[e].visits.end())
                {
                    (*g_)[e].visits.erase(pos);
                }
            }
            else
            {
                (*g_)[e].residual_capacity = 0;
            }
        }
    }
    else
    {
        std::set<unsigned int> uq; // placeholder
        median = reduce_flow(blockedPath, uq, cc, true);
    }
    return median;
}

// calculates the flows and corresponding paths through the graph
void UnitigGraph::assemble(std::string fname, float error_rate)
{
    unsigned int i = 0;
    std::cerr << "Cleaning graph" << std::endl;
    for (unsigned int cc = 0; cc < graphs_.size(); cc++)
    {
        UGraph* g_ = graphs_.at(cc);
        std::string filename = fname + "CC" + std::to_string(cc) + "Graph.dot";
        std::ofstream outfile (filename);
        //contractPaths(cc); // TODO only for debug reasons
        printGraph(outfile, cc);
        cleanGraph(cc, error_rate);
        g_ = graphs_.at(cc);
        std::cerr << "Graph " << cc << ": " << boost::num_vertices(*g_) << " vertices remaining" << std::endl;
        std::cerr << "Calculating paths" << std::endl;
        auto paths = find_paths(cc);
        auto seeds = paths.first;
        auto all_paths = paths.second;
        UEdge seed;
        for (auto e : boost::edges(*g_))
        {
            float path_val = 0.f;
            for (auto&& v : (*g_)[e].visits)
            {
                path_val += all_paths[v - 1];
            }
            if (path_val == 0)
            {
                //TODO
            }
            (*g_)[e].last_visit = 0; // reset last visit
        }
        std::vector<UEdge> sources = get_sources(cc); //get sources of the graph (indegree = 0)
        auto edge_compare = [&](UEdge e1, UEdge e2){ //sort by biggest capacity
            return (*g_)[e1].capacity > (*g_)[e2].capacity;
        };
        std::sort(sources.begin(), sources.end(), edge_compare); // so that we search the highest source first
        while (hasRelevance(cc))
        {
            unvisit(cc);
            UEdge seed;
            bool set = false;
            for (auto&& s : seeds)
            {
                if (!(*g_)[s].visits.empty())
                {
                    seed = s;
                    set = true;
                    break;
                }
            }
            if (!set)
            {
                seed = get_next_source(cc);
            }
            std::cerr << "Fixing flow and finding fattest path" << std::endl;
            std::string filename = fname + "CC" + std::to_string(cc) + "Graph" + std::to_string(i) + ".dot";
            std::ofstream outfile (filename);
            printGraph(outfile, cc);
            std::vector<UEdge> path = fixFlow(seed, cc);
            std::cerr << "Calculating contig " << i << "..." << std::endl;
            auto contig = calculate_contigs(path, all_paths, cc);
            if (contig.first.size() > 150)
            {
                std::cout << ">Contig_" << i << "_flow_" << contig.second << "_cc_" << cc << std::endl;
                std::cout << contig.first << std::endl;
            }
            else
            {
                std::cerr << "Removed short contig" << std::endl;
                std::cerr << contig.first << std::endl;
            }
            i++;
            std::cerr << "Cleaning graph again..." << std::endl;
            cleanPath(path, seeds, cc);
            removeStableSets(cc);
        }
    }
    std::cerr << "Assembly complete" << std::endl;
}

void UnitigGraph::cleanPath(std::vector<UEdge>& path, std::vector<UEdge>& seeds, unsigned int cc)
{
    UGraph* g_ = graphs_.at(cc);
    std::set<UEdge> toDelete;
    for (auto&& e : path)
    {
        if ((*g_)[e].capacity == 0 or (*g_)[e].capacity < thresholds_[cc])
        {
            toDelete.insert(e);
            seeds.erase(std::remove(seeds.begin(), seeds.end(), e), seeds.end()); // dont use deleted edge as seed
        }
    }
    for (auto&& e : toDelete)
    {
        auto source = boost::source(e, *g_);
        boost::remove_edge(e, *g_);
		unsigned int indegree = boost::in_degree(source, *g_);
		unsigned int outdegree = boost::out_degree(source,*g_);
		if (outdegree == 0 and indegree == 0)
		{
			boost::remove_vertex(source, *g_);
		}
    }
}

void UnitigGraph::printGraph(std::ostream& os, unsigned int cc)
{
    UGraph* g_ = graphs_.at(cc);
    typedef std::map<UVertex, int> IndexMap;
    IndexMap mapIndex;
    boost::associative_property_map<IndexMap> propmapIndex(mapIndex);
    uvertex_iter vi, vi_end;
    int i = 1;
    for (boost::tie(vi,vi_end) = boost::vertices(*g_); vi != vi_end; ++vi)
    {
        boost::put(propmapIndex,*vi,i++);
    }
    for (auto e : boost::edges(*g_))
    {
        (*g_)[e].v.visits = (*g_)[e].visits;
    }
    //boost::write_graphviz(os, *g_, boost::make_label_writer(boost::get(&VertexProperties::index,*g_)), boost::make_label_writer(boost::get(&EdgeProperties::name,*g_)), boost::default_writer(), propmapIndex);
    boost::write_graphviz(os, *g_, boost::make_label_writer(boost::get(&VertexProperties::index,*g_)), boost::make_label_writer(boost::get(&EdgeProperties::cap_info,*g_)), boost::default_writer(), propmapIndex);
    //boost::write_graphviz(os, *g_, boost::make_label_writer(boost::get(&VertexProperties::index,*g_)), boost::make_label_writer(boost::get(&EdgeProperties::capacity,*g_)), boost::default_writer(), propmapIndex);
    //boost::write_graphviz(os, *g_, boost::make_label_writer(boost::get(&VertexProperties::index,*g_)), boost::make_label_writer(boost::get(&EdgeProperties::residual_capacity,*g_)), boost::default_writer(), propmapIndex);
    //boost::write_graphviz(os, *g_, boost::make_label_writer(boost::get(&VertexProperties::index,*g_)), boost::make_label_writer(boost::get(&EdgeProperties::v,*g_)), boost::default_writer(), propmapIndex);
    //boost::write_graphviz(os, *g_, boost::make_label_writer(boost::get(&VertexProperties::index,*g_)), boost::make_label_writer(boost::get(&EdgeProperties::distance,*g_)), boost::default_writer(), propmapIndex);
    //boost::write_graphviz(os, *g_, boost::make_label_writer(boost::get(&VertexProperties::index,*g_)), boost::make_label_writer(boost::get(&EdgeProperties::fatness,*g_)), boost::default_writer(), propmapIndex);
}

void UnitigGraph::debug()
{
    std::string filename = "/home/afritz/Documents/Code/Graph/deBruijn_2.0/build/out/Simulated_viromes/kmer_histograms/Cov2.tsv";
    std::ifstream file(filename);
    std::string line;
    std::vector<std::map<unsigned int, unsigned int>> foo;
    std::map<unsigned int, unsigned int> bar;
    while (std::getline(file, line))
    {
        auto pos1 = line.find_first_of("\t ");
        auto pos2 = line.find_last_of("\t ");
        auto kmer = stoi(line.substr(0,pos1));
        auto val = stoi(line.substr(pos2));
        bar.insert(std::make_pair(kmer, val));
    }
    foo.push_back(bar);
    std::string outp = "/home/afritz/Documents/Code/Graph/deBruijn_2.0/build/out/Simulated_viromes/kmer_histograms/";
    get_thresholds(foo, outp, 0.02); 
	// DEBUG
}
