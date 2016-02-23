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
		if (start == '@') // read name. Next line will be the sequence. If quality starts with @, then the next line will as well
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

std::string deBruijnGraph::dfs(const std::string& source, int apply, bool forward, void (*f)(Vertex&, int), bool (*condition)(const Vertex&, int), bool stop = false)
{
	std::stack<std::string> s;
	s.push(source);
	while (s.size() > 0)
	{
		auto curr = s.top();
		s.pop();
		if (condition(graph_[curr], apply))
			if (stop)
				return curr;
			else
				continue;
		else
		{
			f(graph_[curr],apply);
			std::string next;
			const auto& succ = graph_[curr].get_successors();
			for (const auto& c : succ)
			{
				next = curr.substr(1);
				next.push_back(c);
				s.push(next);
			}
		}
	}
	return source;
}

std::string deBruijnGraph::bfs(const std::string& source, int apply, bool forward, void (*f)(Vertex&, int), bool (*condition)(const Vertex&, int), bool stop = false)
{
	std::queue<std::string> q;
	q.push(source);
	while(q.size() > 0)
	{
		auto curr = q.front();
		q.pop();
		if (condition(graph_[curr],graph_[curr].get_successors().size()))
			if (stop)
				{
					graph_[curr].cc = 0;
					return curr;
				}
			else
				continue;
		else
		{
			f(graph_[curr], apply); // apply f to current vertex and with apply as value (e.g. flow or connected component)
			std::string next;
			const auto& succ = graph_[curr].get_successors();
			for (const auto& c : succ)
			{
				next = curr.substr(1);
				next.push_back(c);
				q.push(next);
			}
			if (!forward) // forward only performs forward bfs in directed graph
			{
				const auto& pred = graph_[curr].get_predecessors();
				for (const auto& c : pred)
				{
					next = curr.substr(0,curr.length() - 1);
					next = c + next;
					q.push(next);
				}
			}
		}
	}
	return source; //TODO?
}

int deBruijnGraph::getSize()
{
	return graph_.size();
}

std::vector<std::string> deBruijnGraph::getSources()
{
	std::vector<std::string> sources;
	for (const auto& p : graph_)
		if (p.second.isSource())
			sources.push_back(p.first);
	return sources;
}

std::vector<std::string> deBruijnGraph::getSinks()
{
	std::vector<std::string> sinks;
	for (const auto& p : graph_)
		if (p.second.isSink())
			sinks.push_back(p.first);
	return sinks;
}

std::string deBruijnGraph::find_next_junction(const std::string& source)
{
	const auto& succ = graph_[source].get_successors();
	return bfs(source, 0 /*do we need this?*/, true, [](Vertex& v, int x){v.cc++;},[](const Vertex& v, int x){return (v.cc == x);}, true);
}
