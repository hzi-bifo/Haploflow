#include "deBruijnGraph.h"

deBruijnGraph::deBruijnGraph(unsigned int k) : k_ (k)
{
}

deBruijnGraph::deBruijnGraph(std::string filename, unsigned int k) : k_ (k)
{
	unsigned int i = 0;
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
			i++;
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

std::string deBruijnGraph::reverse_complement(const std::string& kmer)
{
	std::string rc(kmer);
	auto first = rc.begin();
	auto last = rc.end();
	auto lambda = [](char& c){switch (c){case 'A' : return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T' : return 'A'; default : return 'N';}};
	while (first != last and first != --last)
	{
		auto tmp = *first;
		*first = lambda(*last);
		*last = lambda(tmp);
		++first;
	}
	*first = lambda(*first); // make sure k is uneven!
	return rc;
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

std::string deBruijnGraph::dfs(const std::string& source, bool forward, std::function<void(std::string&)> f, std::function<bool(std::string&)> condition, bool stop = false)
{
	std::stack<std::string> s;
	s.push(source);
	while (s.size() > 0)
	{
		auto curr = s.top();
		s.pop();
		if (condition(curr))
			if (stop)
				return curr;
			else
				continue;
		else
		{
			f(curr);
			std::string next;
			const auto& succ = graph_[curr].get_successors();
			for (const auto& c : succ)
			{
				next = curr.substr(1);
				next.push_back(c.first);
				s.push(next);
			}
		}
	}
	return source;
}

std::pair<std::string,unsigned int> deBruijnGraph::bfs(const std::string& source, bool forward, std::function<void(std::string&)> f, std::function<bool(std::string&)> condition, bool stop = false)
{
	unsigned int depth = 0;
	std::queue<std::string> q;
	std::string last = source;
	q.push(source);
	while(q.size() > 0)
	{
		auto curr = q.front();
		q.pop();
		//if (graph_[curr].cc == 0)
			depth++; // depth of the bfs
		if (condition(curr))
		{
			if (stop)
					return std::make_pair(curr,depth);
			else
				continue;
		}
		else
		{
			f(curr); // apply f to current vertex
			std::string next;
			const auto& succ = graph_[curr].get_successors();
			for (const auto& c : succ)
			{
				next = curr.substr(1);
				next.push_back(c.first);
				q.push(next);
				last = next;
			}
			if (!forward) // forward only performs forward bfs in directed graph
			{
				const auto& pred = graph_[curr].get_predecessors();
				for (const auto& c : pred)
				{
					next = curr.substr(0,curr.length() - 1);
					next = c + next;
					q.push(next);
					last = next;
				}
			}
		}
	}
	return std::make_pair(last,depth);
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

std::string deBruijnGraph::find_next_junction(const std::string* source)
{
	const auto& succ = graph_[*source].get_successors();
	if (succ.size() == 0) //sink
	{
		return *source;
	}
	else if (succ.size() == 1)
		return bfs(*source,  
											true, 
											[&, source](std::string& v){graph_[v].source = source;}, 
											[&, source](std::string& v){return (graph_[v].get_successors().size() != 1 or (graph_[v].get_predecessors().size() > 1 and v != *source));},
											true).first;
	else
	{
		unsigned int x = succ.size() - 1;
		const auto& res = bfs(*source, 
																	true, 
																	[&, source](std::string& v){if (graph_[v].source != source){graph_[v].cc++;} graph_[v].source = source;},
																	[&, x](std::string& v){bool ret = (graph_[v].cc == x); if (ret){graph_[v].cc = 0; graph_[v].source = 0;} return ret;}, 
																	true);
		if (res.first == *source)
		{
			auto succ = graph_[*source].get_successors();
			if (res.second > 1)
			{
				// we didnt immediately return, circle found
				std::string neigh = source->substr(1) + succ[0].first;
			}
		}
		return res.first;
	}
	return *source; // we shouldnt end here
}

std::unordered_map<std::string, std::string> deBruijnGraph::find_all_junctions()
{
	std::vector<std::string> sources = getSources();
	auto s_tmp = getSinks();
	std::unordered_set<std::string> sinks(s_tmp.begin(), s_tmp.end());
	std::unordered_map<std::string, std::string> junctions; // maps vertices to their next junction
	for (const auto& source : sources)
	{
		std::string curr = source;
		std::string next = find_next_junction(&curr);
		while (curr != next and junctions.emplace(curr, next).second)
		{
			curr = next;
			next = find_next_junction(&curr);
		}
	}
	//std::cerr << "Junctions: " << junctions.size() << std::endl;
	return junctions;
}

std::vector<std::pair<std::string,unsigned int> > deBruijnGraph::getSequence(const std::pair<std::string,std::string>& junk)
{
	std::vector<std::pair<std::string,unsigned int> > sequences;
	std::string seq = "";
	std::string curr = junk.first;
	std::string next;
	auto neigh = graph_[curr].get_successors();
	for (const auto& n : neigh)
	{
		unsigned int num = n.second;
		std::cout << num << std::endl;
		std::stack<std::string> s;
		s.push(curr);
		while (s.size() > 0)
		{
			auto curr = s.top();
			s.pop();
			if (curr == junk.second)
			{
				sequences.push_back(std::make_pair(seq,num));
				break;
			}
			else
			{
				std::string next;
				const auto& succ = graph_[curr].get_successors();
				bool found_next = false;
				for (const auto& c : succ)
				{
					next = curr.substr(1);
					next.push_back(c.first);
					if (-2 + num <= c.second and 2 + num >= c.second)
					{
						if (found_next)
						{
							std::cout << "double path" << std::endl;
							break;
							// found breakpoint, TODO
						}
						else
						{
							found_next = true;
							seq += c.first;
							num = c.second; //
							s.push(next);
						}
					}
				}
				if (!found_next)
				{
					std::cout << "no path" << std::endl;
					// no suitable neighbour found -> maybe previous part was shared by two strains TODO
				}
			}
		}
	}
	return sequences;
}

std::vector<std::pair<std::string, unsigned int> > deBruijnGraph::getScaffolds(std::unordered_map<std::string, std::string>& junk)
{
	std::vector<std::string> sources = getSources();
	std::vector<std::pair<std::string, unsigned int> > scaffolds;
	for (const auto& source : sources)
	{
		std::string curr = source; 
		std::string next = junk[curr];
		while (graph_[next].get_successors().size() == 0 or junk.find(next) != junk.end())
		{
			auto pair = std::make_pair(curr,next);
			auto seqs = getSequence(pair);
			for (auto& seq : seqs)
			{
				if (seq.second >= 2) //TODO: differentiate between low abundant strain and possible erroneous path 
				{
					scaffolds.push_back(seq);
					//TODO
				}
			}
			if (graph_[next].get_successors().size() == 0)
				break;
			curr = next;
			next = junk[next];
		}
	}
	return scaffolds;
}

void deBruijnGraph::debug()
{
}
