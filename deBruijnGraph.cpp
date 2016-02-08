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

void deBruijnGraph::add_as_neighbour(const std::string& kmer, const char& letter, bool succ)
{
	int i = (succ ? 0 : 4); // successors are in the first 4 slots, predecessors in the second 4
	switch(letter)
	{
		case 'A': graph_[kmer][i] += 1; break;
		case 'C': graph_[kmer][i + 1] += 1; break;
		case 'G': graph_[kmer][i + 2] += 1; break;
		case 'T': graph_[kmer][i + 3] += 1; break;
		default: break; //TODO
	}
}

void deBruijnGraph::split_read(const std::string& line)
{
	std::array<unsigned int,8> init_array = {0,0,0,0,0,0,0,0};
	// the first kmer does not have predecessors, init manually
	std::string kmer = line.substr(0,k_);
	graph_.emplace(kmer,init_array);
	add_as_neighbour(kmer, line[k_],true); // add the k+1st letter as neighbour
	for (unsigned int i = k_ + 1; i < line.length(); i++)
	{
		kmer = line.substr(i - k_,k_); // extract kmer
		graph_.emplace(kmer,init_array); // if not in list, add kmer
		auto next_letter = line[i]; // depending on the next char add an "edge count"
		add_as_neighbour(kmer, next_letter, true);
		auto prev_letter = line[i - k_ - 1];
		add_as_neighbour(kmer, prev_letter, false);
	}
	// this for-loop does not add the final kmer of the read, add manually:
	kmer = line.substr(line.length() - k_, k_);
	graph_.emplace(kmer,init_array); //the last node does not have neighbours, if it already is in the graph, then nothing will change
	add_as_neighbour(kmer, line[line.length() - k_ - 1], false);
}
