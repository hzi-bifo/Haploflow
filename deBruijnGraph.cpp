#include "deBruijnGraph.h"

deBruijnGraph::deBruijnGraph(unsigned int k) : k_ (k)
{
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
	std::cout << graph_.size() << std::endl;
	for (const auto& v: graph_)
	{
		std::cout << v.first << std::endl;
		v.second.print();
		std::cout << std::endl;
	}
}

void deBruijnGraph::split_read(const std::string& line)
{
	// the first kmer does not have predecessors, init manually
	std::string kmer = line.substr(0,k_);
	graph_.emplace(std::pair<std::string, Vertex>(kmer, {})); // add "empty" vertex
	graph_[kmer].add_successor(line[k_]); // add the k+1st letter as neighbour
	for (unsigned int i = k_ + 1; i < line.length(); i++)
	{
		kmer = line.substr(i - k_,k_); // extract kmer
		graph_.emplace(std::pair<std::string, Vertex>(kmer, {})); // if not in list, add kmer
		graph_[kmer].add_successor(line[i]);
		graph_[kmer].add_predecessor(line[i - k_ - 1]);
	}
	// this for-loop does not add the final kmer of the read, add manually:
	kmer = line.substr(line.length() - k_, k_);
	graph_.emplace(std::pair<std::string, Vertex>(kmer, {})); //the last node does not have neighbours, if it already is in the graph, then nothing will change
	graph_[kmer].add_predecessor(line[line.length() - k_ - 1]);
}

std::vector<std::string> deBruijnGraph::get_terminals(bool sink = false)
{
	std::vector<std::string> terminals;
	for (const auto& v : graph_)
	{
		auto& neigh = v.second;
		if (neigh.get_capacity() == 0)
			terminals.push_back(v.first);
	}
	return terminals;
}

void deBruijnGraph::mark_ccs(const std::string& source, unsigned int state)
{
	std::queue<std::string> q;
	q.push(source);
	while(q.size() > 0)
	{
		auto curr = q.front();
		q.pop();
		if (graph_[curr].cc != 0)
			continue;
		else
		{
			graph_[curr].cc = state;
			std::string next;
			const auto& succ = graph_[curr].get_successors();
			for (const auto& c : succ)
			{
				next = curr.substr(1);
				next.push_back(c);
				if (graph_[next].cc == 0)
					q.push(next);
			}
			const auto& pred = graph_[curr].get_predecessors();
			for (const auto& c : pred)
			{
				next = curr.substr(0,curr.length() - 1);
				next = c + next;
				if (graph_[next].cc == 0)
					q.push(next);
			}
		}
	}
}
