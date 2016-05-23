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
		v.print();
		std::cout << std::endl;
	}
}

void deBruijnGraph::split_read(const std::string& line)
{
	// the first kmer does not have predecessors, init manually
	std::string kmer = line.substr(0,k_);
	Vertex toAdd(kmer);
	auto&& v = graph_.emplace(toAdd); // add "empty" vertex
	if (!v.second and v.first->isRC(kmer)) // vertex has been added and was a reverse complement
		v.first->add_predecessor(complement(line[k_])); // if RC(A)->X, then X->A
	else
		v.first->add_successor(line[k_]); // add the k+1st letter as neighbour
	
	for (unsigned int i = k_ + 1; i < line.length(); i++)
	{
		kmer = line.substr(i - k_,k_); // extract kmer
		toAdd = Vertex(kmer);
		v = graph_.emplace(toAdd); // if not in list, add kmer
		if (!v.second and v.first->isRC(kmer))
		{
			v.first->add_predecessor(complement(line[i]));
			v.first->add_successor(complement(line[i - k_ - 1]));
		}
		else
		{
			v.first->add_successor(line[i]);
			v.first->add_predecessor(line[i - k_ - 1]);
		}
	}
	// this for-loop does not add the final kmer of the read, add manually:
	kmer = line.substr(line.length() - k_, k_);
	toAdd = Vertex(kmer);
	v = graph_.emplace(toAdd); //the last node does not have neighbours, if it already is in the graph, then nothing will change
	if (!v.second and v.first->isRC(kmer))
		v.first->add_successor(complement(line[line.length() - k_ - 1]));
	else
		v.first->add_predecessor(line[line.length() - k_ - 1]);
}

std::pair<std::string,unsigned int> deBruijnGraph::bfs(const std::string& source, int* x, std::function<void(const std::string&,int*)> f, std::function<bool(const std::string&, int*)> condition, bool stop = false)
{
	unsigned int depth = 0;
	std::queue<std::string> q;
	std::string last = source;
	q.push(source);
	while(q.size() > 0)
	{
		std::string curr = q.front();
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
			auto&& v = graph_.find(curr);
			auto&& succ = v->get_successors(v->isRC(curr));
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
		if (p.isSource(false) /*or p.isSource(true)*/)
			sources.push_back(p.kmer);
	return sources;
}

std::vector<std::string> deBruijnGraph::getSinks()
{
	std::vector<std::string> sinks;
	for (const auto& p : graph_)
		if (p.isSink(false) /*or p.isSink(true)*/)
			sinks.push_back(p.kmer);
	return sinks;
}

std::string deBruijnGraph::find_next_junction(const std::string* source)
{
	auto v = graph_.find(*source);
	auto&& succ = v->get_successors(v->isRC(*source));
	v->visited = true;
	if (succ.size() == 0) //sink
	{
		return "";
	}
	else if (succ.size() == 1)
	{
		const std::pair<std::string,unsigned int>& res = bfs(*source, 0,
											[&, source](const std::string& v, int* i){
													auto&& s = graph_.find(v);
													s->visited = true;
													}, 
											[&, source](const std::string& v, int* i){
													auto&& s = graph_.find(v);
													return (s->get_successors(s->isRC(v)).size() != 1 or (s->get_predecessors(s->isRC(v)).size() > 1 and v != *source));
													},
											true);
		return res.first;
	}
	else
	{
		int x = 1;
		const std::pair<std::string,unsigned int>& res = bfs(*source, &x,
											[&, source](const std::string& v, int* i){
													auto&& s = graph_.find(v);
													s->visited = true;
													int tmp = s->get_successors(s->isRC(v)).size(); 
													if (tmp > 1) // outdegree > 1: bubble begin
														(*i)++;
													if (tmp == 0)
														(*i)--; // outdegree = 0: sink
													tmp = s->get_predecessors(s->isRC(v)).size();
													if (tmp > 1) // indegree > 1: bubble end
														(*i)--;
												}, 
											[&, source](const std::string& v, int* i){
													auto&& s = graph_.find(v);
													if (s == graph_.end())
														std::cout << v << std::endl;
													return (s->get_successors(s->isRC(v)).size() == 0 or !*i);
												},
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
	auto&& w = graph_.find(sink);
	bool rc = false;
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
			auto&& v = graph_.find(curr);
			for (const auto& n : v->get_successors(v->isRC(curr)))
			{
				std::string next = curr.substr(1);
				next.push_back(n);
				v = graph_.find(next);
				rc = v->isRC(next);
				if (pred.find(next) == pred.end() and next != source and v->capacity(rc) > v->flow(rc))
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
		rc = w->isRC(sink);
		unsigned int max_flow = w->capacity(rc) - w->flow(rc) + 1;
		std::string next = sink;
		while (next != source)
		{
			w = graph_.find(next);
			rc = w->isRC(next);
			max_flow = std::min(max_flow,w->capacity(rc) - w->flow(rc));
			next = pred[next] + next.substr(0,next.size() - 1);
		}
		next = sink;
		while (next != source)
		{
			w = graph_.find(next);
			rc = w->isRC(next);
			path.push_back(pred[next]);
			w->add_flow(rc, max_flow);
			//w->add_flow(!rc, max_flow); //?
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
		std::cout << seq.first << std::endl;
	}
	std::cerr << (clock() - t)/1000000. << std::endl;
}
