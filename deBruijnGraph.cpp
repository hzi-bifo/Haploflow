#include "deBruijnGraph.h"

deBruijnGraph::deBruijnGraph(unsigned int k) : k_ (k)
{
}

deBruijnGraph::deBruijnGraph(std::string filename, bool fasta, unsigned int k) : k_ (k)
{
	unsigned int i = 0;
	// create dBg from FASTA/Q file. FASTA assumed to have only header + sequence in body
	if (!fasta)
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
				i++;
			}
		}
	}
	else
	{
		split_fasta(filename);
	}
}

//currently only for adding fasta sequences
void deBruijnGraph::add_sequence(std::string filename)
{
	split_fasta(filename);
}

void deBruijnGraph::split_fasta(std::string filename)
{
	unsigned int rep_k = 0;
	std::ifstream infile(filename);
	std::string line;
	// first line is header, ignore(?)
	std::string header;
	std::getline(infile, header);
	// save previous read to include the first k-1 signs in the next line if necessary
	std::string prev;
	std::getline(infile, prev);
	if (prev.length() > k_)
		rep_k += split_read(prev);
	else
		std::cerr << "Make sure that linesize it at least k" << std::endl;
	while (std::getline(infile, line))
	{
		std::string to_sep;
		// check how we have to append the previous line to account for line wrapping 
		auto&& linesize = line.length();
		if (linesize > k_ and prev.length())
			to_sep = prev.substr(linesize - k_ + 1) + line;
		else if (linesize)
			to_sep = prev + line;
		else //TODO newlines causes kmer to be added twice, might cause problems
			continue;
		rep_k += split_read(to_sep);
		auto&& fchar = line.front();
		std::string lastk = prev.substr(prev.length() - k_, k_);
		std::string firstk = lastk.substr(1) + fchar;
		auto&& lchar = lastk.front();
		auto&& v = graph_.find(lastk);
		auto&& w = graph_.find(firstk);
		// connect the two lines
		if (v->isRC(lastk))
			v->add_predecessor(complement(fchar));
		else
			v->add_successor(fchar);
		if (w->isRC(firstk))
			w->add_successor(complement(lchar));
		else
			w->add_predecessor(lchar);
		if (v->get_successors(false).size() > 1 or v->get_predecessors(false).size() > 1)
			junctions_.emplace(*v);
		if (w->get_successors(false).size() > 1 or w->get_predecessors(false).size() > 1)
			junctions_.emplace(*w);
		prev = line;	
	}
}

void deBruijnGraph::printGraph()
{
	std::cout << graph_.size() << std::endl;
	for (const auto& v: graph_)
	{
		v.print(false);
		std::cout << std::endl;
	}
}

unsigned int deBruijnGraph::split_read(const std::string& line)
{
	// DEBUG
	unsigned int rep_k = 0;
	// the first kmer does not have predecessors, init manually
	std::string kmer = line.substr(0,k_);
	Vertex toAdd(kmer);
	auto&& v = graph_.emplace(toAdd); // add "empty" vertex
	if (!v.second)//DEBUG
		rep_k++;
	if (!v.second and v.first->isRC(kmer)) // vertex has been added and was a reverse complement
		v.first->add_predecessor(complement(line[k_])); // if RC(A)->X, then X->A
	else
		v.first->add_successor(line[k_]); // add the k+1st letter as neighbour
	if (v.first->get_predecessors(false).size() > 1 or v.first->get_successors(false).size() > 1)
		junctions_.emplace(*v.first);

	for (unsigned int i = k_ + 1; i < line.length(); i++)
	{
		kmer = line.substr(i - k_,k_); // extract kmer
		toAdd = Vertex(kmer);
		v = graph_.emplace(toAdd); // if not in list, add kmer
		if (!v.second)//DEBUG
			rep_k++;
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
		if (v.first->get_predecessors(false).size() > 1 or v.first->get_successors(false).size() > 1)
			junctions_.emplace(*v.first);
	}
	// this for-loop does not add the final kmer of the read, add manually:
	kmer = line.substr(line.length() - k_, k_);
	toAdd = Vertex(kmer);
	v = graph_.emplace(toAdd); //the last node does not have neighbours, if it already is in the graph, then nothing will change
	if (!v.second)//DEBUG
		rep_k++;
	if (!v.second and v.first->isRC(kmer))
		v.first->add_successor(complement(line[line.length() - k_ - 1]));
	else
		v.first->add_predecessor(line[line.length() - k_ - 1]);
	if (v.first->get_predecessors(false).size() > 1 or v.first->get_successors(false).size() > 1)
		junctions_.emplace(*v.first);
	
	return rep_k; //DEBUG
}

template<typename T>
std::pair<std::string,unsigned int> deBruijnGraph::bfs(const std::string& source, T* x, std::function<void(const std::string&,T*)> f, std::function<bool(const std::string&, T*)> condition, bool stop = false)
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

std::vector<std::pair<std::string, unsigned int> > deBruijnGraph::getSequences(const std::string& source, const std::string& sink)
{
	std::vector<std::pair<std::string, unsigned int> > paths;
	unsigned int flow = 0;
	auto&& w = graph_.find(sink);
	bool rc = false;
	// implementation of Edmonds-Karp
	while (true)
	{
		// queue to bfs from
		std::queue<std::string> q;
		q.push(source);
		std::string path = "";
		// store backward path
		std::unordered_map<std::string,char> pred;
		while (q.size() > 0)
		{
			std::string curr = q.front();
			q.pop();
			auto&& v = graph_.find(curr);
			auto&& succ = v->get_successors(v->isRC(curr));
			//create path when there is capacity and we havent been here before
			for (const auto& n : succ)
			{
				std::string next = curr.substr(1);
				next.push_back(n);
				v = graph_.find(next);
				rc = v->isRC(next);
				// if next in pred: we have found a cycle before finding sink
				if (pred.find(next) == pred.end() and v->capacity(rc) > v->flow(rc))
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
		unsigned int max_flow = w->capacity(rc) + 1;
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
			next = pred[next] + next.substr(0,next.size() - 1);
		}
		std::reverse(path.begin(),path.end());
		//path += sink; // careful if "glueing" is about to take place!
		paths.push_back(std::make_pair(path,max_flow));
		flow += max_flow;
	}
	return paths;
}

std::string deBruijnGraph::find_next_junction(const std::string& source)
{
	bfs<unsigned int>(source,0,
					[&,source](const std::string& source, unsigned int* i)
					{
						auto&& v = graph_.find(source);
						v->visited = true;
					},
					[&,source](const std::string& source, unsigned int* i)
					{
						auto&& v = graph_.find(source);
						if (v.get_successors(v->isRC(source)).size() > 1 or v.get_predecessors(v->isRC(source)).size() > 1)
						{
							//junction found

						}
						return v->visited;
					},
					true);
	return "LUL";
}

std::string deBruijnGraph::make_graph()
{
	std::string ret;
	return ret;
}

void shrink(const std::string source)
{
		
}


void deBruijnGraph::debug()
{
	std::cerr << "Vertices: " << getSize() << std::endl;
	clock_t t = clock();
	std::vector<std::string> sources = getSources();
	std::vector<std::string> sinks = getSinks();
	std::cerr << sources.size() << " sources found" << std::endl;
	std::cerr << sinks.size() << " sinks found" << std::endl;
	auto&& path = getSequences(sources[0],sinks[0]);
	for (const auto& p : path)
		std::cout << p.first << std::endl;
	std::cerr << (clock() - t)/1000000. << std::endl;
}
