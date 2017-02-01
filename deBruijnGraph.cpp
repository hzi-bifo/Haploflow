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

// returns Sequence in graph, throws exception if not in graph
const Sequence& deBruijnGraph::getSequence(const std::string& kmer)
{
	Sequence seq(kmer);
	auto&& ret = graph_.find(seq);
	if (ret != graph_.end())
		return ret->first;
	else
	{
		std::cerr << "kmer not in graph: " << kmer << std::endl;
		throw;
	}
}

Vertex* deBruijnGraph::getVertex(const std::string& kmer)
{
	if (kmer.length() != k_)
		return 0;
	else
	{
		Sequence seq(kmer);
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
