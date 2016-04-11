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

std::pair<std::string,unsigned int> deBruijnGraph::bfs(const std::string& source, int* x, std::function<void(std::string&,int*)> f, std::function<bool(std::string&, int*)> condition, bool stop = false)
{
	unsigned int depth = 0;
	std::queue<std::string> q;
	std::string last = source;
	q.push(source);
	while(q.size() > 0)
	{
		auto curr = q.front();
		q.pop();
		depth++; // depth of the bfs
		if (condition(curr,x))
		{
			if (stop)
					return std::make_pair(curr,depth);
			else
				continue;
		}
		else
		{
			f(curr,x); // apply f to current vertex
			std::string next;
			const auto& succ = graph_[curr].get_successors();
			for (const auto& c : succ)
			{
				next = curr.substr(1);
				next.push_back(c);
				q.push(next);
				last = next;
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
	graph_[*source].visited = true;
	if (succ.size() == 0) //sink
	{
		return "";
	}
	else if (succ.size() == 1)
	{
		const auto& res = bfs(*source, 0,
											[&, source](std::string& v, int* i){graph_[v].visited = true;}, 
											[&, source](std::string& v, int* i){return (graph_[v].get_successors().size() != 1 or (graph_[v].get_predecessors().size() > 1 and v != *source));},
											true);
		return res.first;
	}
	else
	{
		int x = 1;
		const auto& res = bfs(*source, &x,
											[&, source](std::string& v, int* i){
																						graph_[v].visited = true;
																						int tmp = graph_[v].get_successors().size(); 
																						if (tmp > 1) // outdegree > 1: bubble begin
																							(*i)++;
																						if (tmp == 0)
																							(*i)--; // outdegree = 0: sink
																						tmp = graph_[v].get_predecessors().size();
																						if (tmp > 1) // indegree > 1: bubble end
																							(*i)--;}, 
											[&, source](std::string& v, int* i){return (graph_[v].get_successors().size() == 0 or !*i);},
											true);
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
		while (next != "" and junctions.emplace(curr, next).second)
		{
			curr = next;
			next = find_next_junction(&curr);
		}
	}
	//std::cerr << "Junctions: " << junctions.size() << std::endl;
	return junctions;
}

std::vector<std::pair<std::string, unsigned int> > deBruijnGraph::getSequences (const std::string& source, const std::string& sink)
{
	std::vector<std::pair<std::string, unsigned int> > paths;
	unsigned int flow = 0;
	while (true)
	{
		std::queue<std::string> q;
		q.push(source);
		std::string path = "";
		std::unordered_map<std::string,char> pred;
		while (q.size() > 0)
		{
			std::string curr = q.front();
			q.pop();
			for (const auto& n : graph_[curr].get_successors())
			{
				std::string next = curr.substr(1);
				next.push_back(n);
				if (pred.find(next) == pred.end() and next != source and graph_[next].capacity() > graph_[next].flow)
				{
					pred[next] = curr[0];
					q.push(next);
				}
			}
		}
		if (pred.find(sink) == pred.end())
		{
			break;
		}
		unsigned int max_flow = graph_[sink].capacity() - graph_[sink].flow + 1;
		std::string next = sink;
		while (next != source)
		{
			max_flow = std::min(max_flow,graph_[next].capacity() - graph_[next].flow);
			next = pred[next] + next.substr(0,next.size() - 1);
		}
		next = sink;
		while (next != source)
		{
			path.push_back(pred[next]);
			graph_[next].flow += max_flow;
			next = pred[next] + next.substr(0,next.size() - 1);
		}
		std::reverse(path.begin(),path.end());
		path = path + sink;
		paths.push_back(std::make_pair(path,max_flow));
		flow += max_flow;
	}
	return paths;
}

void deBruijnGraph::debug()
{
	std::cerr << getSize() << std::endl;
		
	clock_t t = clock();
	auto c = find_all_junctions();
	std::cerr << (clock() - t)/1000000. << std::endl;
	t = clock();
	std::vector<std::pair<std::string,unsigned int> > sequences;
	for (const auto& p : c)
	{
		auto s = getSequences(p.first, p.second);
		sequences.insert(sequences.end(),s.begin(), s.end());
	}
	unsigned int i = 0;
	for (const auto& seq : sequences)
	{
		std::cout << "> Contig_" << i++ << std::endl;
		std::cout << seq.first << " (cov " << seq.second << ")" << std::endl;
	}
	std::cerr << (clock() - t)/1000000. << std::endl;
}
