#include "deBruijnGraph.h"

deBruijnGraph::deBruijnGraph(unsigned int k) : k_ (k)
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
			if (line.length() >= k_)
			{
				split_read(line);
				i++;
			}
			//else 
			//	std::cerr << "Read length less than k, skipping" << std::endl; 
		}
	}
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
				v.scc = index;
				index++;
				v.onStack = true;
				visit_stack.push(*cseq);
			}

			if (child < children)
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
					recursion_stack.push(std::make_pair(std::make_pair(curr,next),0));
				}
				else if (w.onStack)
				{
					v.scc = std::min(v.scc, w.index);
				}
			}
			else  // all children have been visited in the stack
			{
				if (prev != "") // is not the first searched
				{
					const Sequence* pseq = getSequence(prev);
					Vertex& w = graph_[*pseq];
					w.scc = std::min(w.scc, v.scc);
					if (w.index == w.scc)
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
	for (auto&& v : graph_)
	{
		if (v.first == "AGATAGAGTGATGGTATCACCTTTGGCTGTGACATGGTGGA"
		or v.first == "TCCACCATGTCACAGCCAAAGGTGATACCATCACTCTATCT")
		{
			v.second.print(true);
		}
	}
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
