#include "deBruijnGraph.h"

deBruijnGraph::deBruijnGraph(unsigned int k)
{
	k_ = k;
}

deBruijnGraph::deBruijnGraph(std::string filename, unsigned int k) : k_ (k)
{
	std::ifstream infile(filename);
	std::string line;
	bool next_read = false;
	while (std::getline(infile,line))
	{
		const auto& start = line.front();
		if (start == '@') // read name. Next line will be the sequence
		{
			next_read = true;
		}
		else if (next_read)
		{
			next_read = false;
			split_read(line);
		}
	}
}

void deBruijnGraph::printGraph()
{
	for (const auto& v: graph_)
	{
		std::cout << v.first << ": ";
		for (const auto& i: v.second)
		{
			std::cout << i << " ";
		}
		std::cout << std::endl;
	}
}

void deBruijnGraph::split_read(std::string line)
{
	std::array<unsigned int,4> init_array = {0,0,0,0};
	for (unsigned int i = k_; i < line.length(); i++)
	{
		std::string kmer = line.substr(i - k_,k_); // extract kmer
		graph_.emplace(kmer,init_array); // if not in list, add kmer
		switch (line[i]) // depending on the next char add an "edge count"
		{
			case 'A': graph_[kmer][0] += 1; break;
			case 'C': graph_[kmer][1] += 1; break;
			case 'G': graph_[kmer][2] += 1; break;
			case 'T': graph_[kmer][3] += 1; break;
			default: break; //TODO
		}
	}
	// this for-loop does not add the final kmer of the read, add manually:
	graph_.emplace(line.substr(line.length() - k_, k_),init_array); //the last node does not have neighbours, if it already is in the graph, then nothing will change
}
