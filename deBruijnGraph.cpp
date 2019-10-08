#include "deBruijnGraph.h"

deBruijnGraph::deBruijnGraph(std::string filename)
{
    std::ifstream graph_file(filename);
    short counter = -2;
    std::string line;
        
    std::string sequence;
    int a_in; int c_in; int g_in; int t_in;
    int a_out; int c_out; int g_out; int t_out;
    unsigned int starts_with; unsigned int ends_with;

    while (std::getline(graph_file, line))
    {    
        if (counter < 0) // read header
        {
            std::stringstream ss(line);
            while(std::getline(ss, line,'\t'))
            {
                if (counter == -2)
                    k_ = stoi(line);
                else if (counter == -1)
                    read_length_ = stoi(line);
                counter++;
            }
            continue;
        }
        if (line.empty()) // empty lines between vertices
            continue;
        std::stringstream ss(line);
        while (std::getline(ss, line, '\t')) // split the line at tabs
        {
            switch(counter){    
                case 0: sequence = line; break;
                case 1: a_in = stoi(line); break;
                case 2: c_in = stoi(line); break;
                case 3: g_in = stoi(line); break;
                case 4: t_in = stoi(line); break;
                case 5: a_out = stoi(line); break;
                case 6: c_out = stoi(line); break;
                case 7: g_out = stoi(line); break;
                case 8: t_out = stoi(line); break;
                case 9: break;
                case 10: break;
                case 11: starts_with = stoi(line); break; // these don't set anything
                case 12: ends_with = stoi(line); break;
                default: counter = 0; break; 
            }
            counter++; counter %= 13; 
        }
        if (!counter)
        {
            Sequence s(sequence);
            Vertex v(a_in, c_in, g_in, t_in, a_out, c_out, g_out, t_out, std::make_pair(starts_with, ends_with));
            graph_.emplace(s,v);
        }
    }
}

deBruijnGraph::deBruijnGraph(unsigned int k, std::unordered_map<Sequence, Vertex> vertices) : graph_(vertices), k_ (k)
{
}

deBruijnGraph::deBruijnGraph(std::string filename, unsigned int k) : k_ (k)
{
	unsigned int i = 0;
	// create dBg from FASTA/Q file. Currently expects one-lined sequences 
	std::ifstream infile(filename);
	std::string line;
	bool next_read = false;
	while (std::getline(infile,line))
	{
		const auto& start = line.front();
		if (start == '@' or start == '>') // read name. Next line will be the sequence. If quality starts with @, then the next line will as well
		{
			next_read = true;
		}
		else if (next_read)
		{
			next_read = false;
			if (line.length() >= k_ and line.find_first_not_of("ACGT") == std::string::npos)
			{
				split_read(line);
				i++;
			}
			//else 
			//	std::cerr << "Read length less than k, skipping" << std::endl; 
		}
	}
}

std::ostream& operator<<(std::ostream& os, const deBruijnGraph& dbg)
{
    os << dbg.k_ << '\t' << dbg.read_length_ << std::endl;
    for (const auto& v : dbg.graph_)
    {
        os << v.first << std::endl;
        os << v.second << std::endl;
    }
    return os;
}

void deBruijnGraph::printGraph() const
{
	std::cout << graph_.size() << std::endl;
	for (const auto& v: graph_)
	{	
		std::cout << v.first.get_kmer() << std::endl;
		v.second.print(false);
		std::cout << std::endl;
	}
}

unsigned int deBruijnGraph::getK() const
{
	return k_;
}

unsigned int deBruijnGraph::split_read(const std::string& line)
{
	// the first kmer does not have predecessors, init manually
	std::string kmer = line.substr(0,k_);
	Sequence toAdd(kmer);
	auto&& v = graph_.emplace(toAdd,Vertex());
	if (!v.second and v.first->first != kmer) // vertex has been added and was a reverse complement
	{
		v.first->second.add_predecessor(complement(line[k_])); // if RC(A)->X, then X->A
	}
	else
	{
		v.first->second.add_successor(line[k_]); // add the k+1st letter as neighbour
	}
	v.first->second.read_start(); // the read started with vertex v

	for (unsigned int i = k_ + 1; i < line.length(); i++)
	{
		kmer = line.substr(i - k_,k_); // extract kmer
		toAdd = Sequence(kmer);
		v = graph_.emplace(toAdd,Vertex()); // if not in list, add kmer
		if (!v.second and v.first->first != kmer)
		{
			v.first->second.add_predecessor(complement(line[i]));
			v.first->second.add_successor(complement(line[i - k_ - 1]));
		}
		else
		{
			v.first->second.add_successor(line[i]);
			v.first->second.add_predecessor(line[i - k_ - 1]);
		}
	}
	// this for-loop does not add the final kmer of the read, add manually:
	kmer = line.substr(line.length() - k_, k_);
	toAdd = Sequence(kmer);
	v = graph_.emplace(toAdd,Vertex()); //the last node does not have neighbours, if it already is in the graph, then nothing will change
	if (!v.second and v.first->first != kmer)
	{
		v.first->second.add_successor(complement(line[line.length() - k_ - 1]));
	}
	else
	{
		v.first->second.add_predecessor(line[line.length() - k_ - 1]);
	}
    v.first->second.read_end(); // the read ended with this vertex
	return 0;
}

int deBruijnGraph::getSize() const
{
	return graph_.size();
}

void deBruijnGraph::markCycles() //non-recusrive tarjan implementation
{
	std::stack<std::pair<std::pair<std::string, std::string>, unsigned int> > recursion_stack;
	std::stack<Sequence> visit_stack;
	unsigned int index = 1;
	for (auto& p : graph_)
	{
		Sequence s = p.first;
		Vertex& v = p.second;
		if (v.index == 0) // scc hasnt been set
		{
			recursion_stack.push(std::make_pair(std::make_pair("",s.get_kmer()),0));
		}
		while (recursion_stack.size() > 0)
		{
			auto next_element = recursion_stack.top();
			recursion_stack.pop();
			
			std::string prev = next_element.first.first; // where we came from
			std::string curr = next_element.first.second;
			const Sequence* cseq = getSequence(curr);
			Vertex& v = graph_[*cseq];
			bool reverse = (*cseq != curr);
			
			unsigned int child = next_element.second;
			unsigned int children = reverse ? v.get_predecessors().size() : v.get_successors().size();
			
			if (!child) // this vertex is visited the first time this run
			{
				v.index = index;
				v.cc = index;
				index++;
				v.onStack = true;
				visit_stack.push(*cseq);
			}

			if (child < children) // still have to search at leat one child
			{
				std::string next("");
				const Sequence* succ = nullptr;
				
				if (!reverse)
				{
					auto successors = v.get_successors();
					next = curr.substr(1) + successors[child];
					succ = getSequence(next);
				}
				else
				{
					auto predecessors = v.get_predecessors();
					next = curr.substr(1) + complement(predecessors[child]); // TODO
					succ = getSequence(next);
				}
				Vertex& w = graph_[*succ];

				recursion_stack.push(std::make_pair(std::make_pair(prev,curr),++child)); // visit the next child
				if (w.index == 0)
				{
					recursion_stack.push(std::make_pair(std::make_pair(curr,next),0)); // "recurse" on current child
				}
				else if (w.onStack)
				{
					v.cc = std::min(v.cc, w.index);
				}
			}
			else  // all children have been visited in the stack
			{
				if (prev != "") // is not the first searched
				{
					const Sequence* pseq = getSequence(prev);
					Vertex& w = graph_[*pseq];
					w.cc = std::min(w.cc, v.cc);
					if (w.index == w.cc)
					{
						Sequence scc = visit_stack.top();
						do
						{
							visit_stack.pop();
							Vertex& path = graph_[scc];
							path.onStack = false;
							scc = visit_stack.top();
						} while (scc != *pseq);
					}
				}
				else
				{
					Sequence top = visit_stack.top();
					Vertex& tv = graph_[top];
					tv.onStack = false;
					visit_stack.pop();
				}
			}
		}
	}
}

unsigned int deBruijnGraph::split_ccs()
{
    unsigned int cc = 1;
    for (auto&& it = graph_.begin(); it != graph_.end(); ++it)
    {
        Sequence s = (*it).first;
        Vertex v = (*it).second;
        if (v.cc == 0)
        {
            auto members = dfs(s, cc++);
        }
    }
    std::cerr << cc << " total connected components" << std::endl;
    return cc;
}
    
std::vector<const Sequence*> deBruijnGraph::dfs(Sequence& s, unsigned int cc)
{
    const std::string kmer = s.get_kmer();
    std::stack<std::pair<const Sequence*, Vertex*>> to_search;
    const Sequence* seq = getSequence(kmer);
    Vertex* v = getVertex(kmer);
    to_search.push(std::make_pair(seq, v));

    std::vector<const Sequence*> members;
    while (!to_search.empty())
    {
        auto curr = to_search.top();
        to_search.pop();
        Vertex* v = curr.second;
        const Sequence* s = curr.first;
        if (v->cc != 0) // has been sarched before
        {
            continue;
        }
        v->cc = cc;
        members.push_back(s);
        std::string seq = s->get_kmer();
        std::string next = "";
        const Sequence* nextS = getSequence(seq); // to check whether sequence is reverse complement or not
        Vertex* nextV;
        bool reverse = (*nextS != seq);

        for (auto& n : v->get_successors())
        {
            if (!reverse)
            {
                next = seq.substr(1) + n;
            }
            else
            {
                next = complement(n) + seq.substr(0,seq.length() - 1);
            }
            nextV = getVertex(next);
            if (nextV == 0)
            {
                std::cerr << (reverse ? "reverse" : "not reverse") << std::endl;
                v->print(true);
                std::cerr << next << " not in graph" << std::endl;
                continue;
            }
            else if (nextV->cc != 0)
            {
                continue;
            }
            nextS = getSequence(next);
            to_search.push(std::make_pair(nextS, nextV));
        }
        for (auto& n : v->get_predecessors())
        {
            if (!reverse)
            {
                next = n + seq.substr(0, seq.length() - 1);
            }
            else
            {
                next = seq.substr(1) + complement(n);
            }
            nextV = getVertex(next);
            if (nextV == 0)
            {
                std::cerr << (reverse ? "reverse" : "not reverse") << std::endl;
                v->print(true);
                std::cerr << next << " not in graph" << std::endl;
                continue;
            }
            else if (nextV->cc != 0)
            {
                continue;
            }
            nextS = getSequence(next);
            to_search.push(std::make_pair(nextS, nextV));
        }
    }
    return members;
}

// calculates some metrics on the de bruijn graph used for estimating cutoffs etc
std::vector<std::map<unsigned int, unsigned int>> deBruijnGraph::coverageDistribution(unsigned int ccs) const
{
    std::vector<std::map<unsigned int, unsigned int>> all_coverages(ccs - 1);
	for (const auto& p : graph_)
	{
		auto& v = p.second;
        auto cc = v.cc - 1; //ccs start at 1
		unsigned int coverage = std::max(v.get_total_in_coverage(), v.get_total_out_coverage());
        //unsigned int coverage = v.get_total_in_coverage() + v.get_total_out_coverage();
		if (all_coverages[cc].find(coverage) != all_coverages[cc].end())
		{
			all_coverages[cc][coverage]++;
		}
		else
		{
			all_coverages[cc][coverage] = 1;
		}
	}
	return all_coverages;
}

std::vector<std::string> deBruijnGraph::getSources() const
{
	std::vector<std::string> sources;
	for (const auto& p : graph_)
		if (p.second.isSource())
			sources.push_back(p.first.get_kmer());
	return sources;
}

std::vector<std::string> deBruijnGraph::getSinks() const
{
	std::vector<std::string> sinks;
	for (const auto& p : graph_)
		if (p.second.isSink())
			sinks.push_back(p.first.get_kmer());
	return sinks;
}

std::pair<std::vector<Sequence>, std::vector<Sequence> > deBruijnGraph::getJunctions() const
{
	std::vector<Sequence> out_unbalanced;
	std::vector<Sequence> in_unbalanced;
	for (auto&& p : graph_)
	{
		unsigned int succ = p.second.get_successors().size();
		unsigned int pred = p.second.get_predecessors().size();
		if (succ > pred)
			out_unbalanced.push_back(p.first);
		else if (pred > succ)
			in_unbalanced.push_back(p.first);
		else if (pred == succ and pred > 1)
		{
			out_unbalanced.push_back(p.first);
			//in_unbalanced.push_back(p.first); TODO oBdA?
		}
	}
	return std::make_pair(out_unbalanced,in_unbalanced);

}

// returns Sequence in graph, returns nullptr if not in graph
const Sequence* deBruijnGraph::getSequence(const std::string& kmer)
{
	Sequence seq(kmer);
	auto&& ret = graph_.find(seq);
	if (ret != graph_.end())
		return &(ret->first);
	else
	{
		return nullptr;
	}
}

Vertex* deBruijnGraph::getVertex(const std::string& kmer)
{
	if (kmer.length() != k_)
		return nullptr;
	else
	{
		Sequence seq(kmer);
		//Sequence seq = getSequence(kmer); //will be the same vertex but constructor is faster
		try
		{
			auto&& v = graph_.at(seq);
			return &v;
		}
		catch (std::out_of_range e)
		{
			return 0;
		}
	}
}

void deBruijnGraph::debug()
{
	//std::cerr << "Vertices: " << getSize() << std::endl;
	/*std::unordered_map<unsigned int, unsigned int> sccs;
	for (auto& p : graph_)
	{
		unsigned int scc = p.second.cc;
		if (sccs.find(scc) != sccs.end())
		{
			sccs[scc]++;
			//std::cout << p.first << " (" << scc << ")" << std::endl;
		}
		else
		{
			sccs[scc] = 2;
		}
	}
	std::cout << sccs.size() << " SCCs" << std::endl;
	unsigned int small_sccs = 0;
	for (const auto& p : sccs)
	{
		if (p.second > 2 )
			std::cout << p.first << " size " << p.second << std::endl;
		else
			small_sccs++;
	}
	std::cout << small_sccs << " SCCs of size 2" << std::endl;*/
	//std::cout << mean << std::endl;
	/*std::cerr << "Vertices: " << getSize() << std::endl;
	clock_t t = clock();
	std::vector<std::string> sources = getSources();
	std::vector<std::string> sinks = getSinks();
	//for (const auto& s : sources)
	//	std::cerr << s << " (source)" << std::endl;
	//for (const auto& t : sinks)
	//	std::cerr << t << " (sink)" << std::endl;
	std::cerr << (clock() - t)/1000000. << std::endl;*/
}
