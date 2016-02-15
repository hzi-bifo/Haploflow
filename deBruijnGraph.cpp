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
	std::cout << graph_.size() << std::endl;
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
		default: graph_[kmer][i] += 1; break; //TODO (currently: add 'A')
	}
}

void deBruijnGraph::split_read(const std::string& line)
{
	std::array<unsigned int,12> init_array = {0,0,0,0,0,0,0,0,0,0,0,0};
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

std::vector<std::string> deBruijnGraph::get_terminals(bool sink = false)
{
	std::vector<std::string> terminals;
	for (const auto& v : graph_)
	{
		auto& neigh = v.second;
		int i = (sink ? 0 : 4); // sinks do not have edges in [0] to [3], sources in [4] to [7]
		if (neigh[i] + neigh[i + 1] + neigh[i + 2] + neigh[i + 3] == 0)
			terminals.push_back(v.first);
	}
	return terminals;
}

void deBruijnGraph::add_back(int i, std::string& next)
{
	if (i == 0)
		next.append("A");
	else if (i == 1)
		next.append("C");
	else if (i == 2)
		next.append("G");
	else if (i == 3)
		next.append("T");
}

void deBruijnGraph::add_front(int i, std::string& next)
{
	if (i == 4)
		next = "A" + next;
	else if (i == 5)
		next = "C" + next;
	else if (i == 6)
		next = "G" + next;
	else if (i == 7)
		next = "T" + next;
}

bool deBruijnGraph::mark_ccs(const std::string& source, unsigned int state, bool forward = true)
{
	std::queue<std::string> q;
	q.push(source);
	while (q.size() > 0)
	{
		auto curr = q.front();
		q.pop();
		if (graph_[curr][8] != 0 and forward)
		{
			mark_ccs(curr, graph_[curr][8], false); // backtracking
			return false;
		}
		else if (graph_[curr][8] == state)
			return false;
		else
			graph_[curr][8] = state;
		std::string next;
		if (forward)
		{
			for (int i = 0; i < 4; i++)
			{
				const auto& n = graph_[curr][i];
				if (n == 0)
					continue;
				else
				{
					next = curr.substr(1);
					add_back(i,next);
				}
				q.push(next);
			}
		}
		else
		{
			for (int i = 4; i < 8; i++)
			{
				const auto& n = graph_[curr][i];
				if (n == 0)
					continue;
				else
				{
					next = curr.substr(0,curr.length() - 1);
					add_front(i,next);
				}
				q.push(next);
			}
		}
	}
	return true;
}

// the node is a source, we search forward until we find a loop or a sink
std::string deBruijnGraph::search(const std::string& source, bool isSource)
{
	std::string seq;
	int m = 0;
	if (isSource)
		seq = source;
	else
	{
		m += 4;
		std::string rev = source;
		std::reverse(rev.begin(), rev.end());
		seq = rev; // so we can append to the back of the string
	}
	std::string next = source;
	bool added = true;
	while(added)
	{
		added = false;
		for (int i = m; i < m+4; i++)
		{
			const auto& n = graph_[next][i];
			if (n == 0)
				continue;
			if (i == 0 or i == 4)
				seq += "A";
			else if (i == 1 or i == 5)
				seq += "C";
			else if (i == 2 or i == 6)
				seq += "G";
			else if (i == 3 or i == 7)
				seq += "T";
			graph_[next][10] = 1; // mark as visited
			next = (isSource ? next.substr(1) : next.substr(0,next.length() - 1));
			isSource ? add_back(i,next) : add_front(i,next);
			added = true;
			break; // just add one possiblity
		}
		if (graph_[next][10])
			break; // we're in a cycle
	}
	if (!isSource)
		std::reverse(seq.begin(),seq.end()); // necessary?
	return seq;
}

std::string deBruijnGraph::extractSequence(const std::string& source, int i = 0)
{
	bool isSource = (graph_[source][4] + graph_[source][5] + graph_[source][6] + graph_[source][7] == 0);
	bool isSink = (graph_[source][0] + graph_[source][1] + graph_[source][2] + graph_[source][3] == 0); // check whether the node is a source or a sink
	if (isSource and isSink)
		return "ERROR";
	else
	{
		std::string out = search(source, !isSink);
		/*debug */
		std::cout << ">contig_" << i << std::endl;
		std::cout << out << std::endl;
		/*debug end*/
		return out;
	}
}

std::string deBruijnGraph::edmonds_karp_single(const std::string& source, const std::string& sink)
{
	std::queue<std::string> q;
	std::string path = source;
	q.push(source);
	while (!q.empty())
	{
		std::string curr = q.front();
		q.pop();
		for (unsigned int i = 0; i < 4; i++)
		{
			if (graph_[curr][i] != 0)
			{
				std::string next = curr.substr(1);
				add_back(i,next);
				if (!graph_[next][9] and next != source and graph_[curr][10] < graph_[curr][11])
				{
					add_back(i,path);
					graph_[next][9] = 1;
					q.push(next);
				}
			}
		}
	}
	if (graph_[sink][9] == 0)
		return; // no path from source to sink
	unsigned int send_flow = std::numeric_limits<unsigned int>::max(); //inf
	
	int j = 0;
	for (std::string s = sink; s != source;) // find the "augmenting path"
	{
		send_flow = std::min(send_flow, graph_[s][11] - graph_[s][10]);
		s = path.substr(path.length() - s.length() - j - 1, s.length());
		j++;
	}
	for (std::string s = sink; s != source;) // let it floooooow
	{
		graph_[s][10] += send_flow;
		s = path.substr(path.length() - s.length() - j - 1, s.length());
		j++;
	}
	return path;
}

void deBruijnGraph::resetVisits()
{
	for (auto& v : graph_)
		v.second[9] = 0;
}
