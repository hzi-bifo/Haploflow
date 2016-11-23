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
	split_read(prev);
	while (std::getline(infile, line))
	{
		std::string to_sep;
		// check how we have to append the previous line to account for line wrapping 
		auto&& linesize = line.length();
		if (linesize > k_ and prev.length())
			to_sep = prev.substr(linesize - k_) + line;
		else if (linesize)
			to_sep = prev + line;
		else
			continue;
		rep_k += split_read(to_sep);
		prev = line;
	}
	//DEBUG
	std::cerr << rep_k << " kmers appeared >1 times" << std::endl;
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
	// the first kmer does not have predecessors, init manually
	unsigned int rep_k = 0;
	std::string kmer = line.substr(0,k_);
	Vertex toAdd(kmer);
	auto&& v = graph_.emplace(toAdd); // add "empty" vertex
	if (!v.second)//DEBUG
		rep_k++;
	if (!v.second and v.first->isRC(kmer)) // vertex has been added and was a reverse complement
		v.first->add_predecessor(complement(line[k_])); // if RC(A)->X, then X->A
	else
		v.first->add_successor(line[k_]); // add the k+1st letter as neighbour

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
	return rep_k; //DEBUG
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
													return (s->get_successors(s->isRC(v)).size() != 1 or (s->get_predecessors(s->isRC(v)).size() != 1 and v != *source));
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
	std::unordered_map<std::string, std::string> junctions; // maps vertices to their next junction
	for (const auto& source : sources)
	{
		std::string curr = source;
		std::string next = find_next_junction(&curr);
		while (next != "" and junctions.find(graph_.find(curr)->rc()) == junctions.end() and junctions.emplace(curr, next).second)
		{
			curr = next;
			next = find_next_junction(&curr);
		}
	}
	//std::cerr << "Junctions: " << junctions.size() << std::endl;
	return junctions;
}

// returns all possible sequences between a single source and sink (maximum 4?)
std::vector<std::pair<std::string, unsigned int> > deBruijnGraph::getSequences (const std::string& source, const std::string& sink)
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

// glues together shorter contigs to "unitigs"
std::pair<std::string, unsigned int> deBruijnGraph::glue(const std::string& source, const std::unordered_map<std::string, std::string>& junctions)
{
	std::string contig = "";
	unsigned int coverage = 0;
	std::string nsource = source; // needs source to start with
	std::string nsink = junctions.at(nsource);
	while (true) // until the sink is a "real" sink (or its RC)
	{
		auto seqs = getSequences(nsource, nsink); //internally uses flow to determine sequence between current source/sink pair
		//std::cout << seqs.size() << std::endl; // number of following sequences
		if (seqs.size() == 1) // single path, glueing is easy
		{
			contig += seqs[0].first;
			if (seqs[0].second < coverage) // this is currently for debugging only, TODO proper coverage handling
				coverage = seqs[0].second;
		}
		else if (seqs.size() > 1) // multiple pathes -> "haplotype" detected -> the magic happens here
		{
			unsigned int max_cov = 0;
			unsigned int i = 0;
			unsigned int pos = 0;
			for (const auto& seq : seqs)
			{
				if (seq.second > max_cov)
				{
					max_cov = seq.second;
					pos = i;
				}
				i++;
			}
			if (seqs[pos].second < coverage)
				coverage = max_cov;
			contig += seqs[pos].first;
		}
		auto v = graph_.find(nsink);
		// stop if sink (only "real sink" = no outgoing)
		if (v->isSink(false))
			break;
		nsource = nsink;
		try
		{
			nsink = junctions.at(nsource);
		}
		// if nsink only has a reverse complement neighbour -> what we find will be the reverse strand
		// nsink is not a "real" sink
		// do not append unless we find diversions on "the way back"?
		catch (const std::out_of_range& err) 
		{
			break; //TODO
			//auto v = graph_.find(nsource);
			//nsink = junctions.at(v->rc());
		}
	}
	return std::make_pair(contig,coverage);
}

std::string deBruijnGraph::make_graph()
{
	std::string ret;
	return ret;
}

void deBruijnGraph::debug()
{
	std::cerr << "Vertices: " << getSize() << std::endl;
		
	clock_t t = clock();
	auto c = find_all_junctions();
	std::cerr << "Partitioned into " << c.size() << " components" << std::endl;
	std::cerr << (clock() - t)/1000000. << std::endl;
	t = clock();
	std::vector<std::string> sources = getSources();
	std::cerr << sources.size() << " sources found" << std::endl;
	// This isnt working as intended
	/*unsigned int i = 0;
	std::vector<std::pair<std::string,unsigned int> > sequences;
	for (const auto& p : sources)
	{
		auto seq = glue(p,c);
		if (seq.first.size() > 0) // otherwise we found a reverse complement
		{
			sequences.push_back(seq);
			std::cout << ">Contig_" << i++ << std::endl;
			std::cout << seq.first << std::endl;
		}
	}
	*/
	/*std::vector<std::pair<std::string,unsigned int> > sequences;
	i = 0;
	for (const auto& p : c)
	{
		auto s = getSequences(p.first, p.second);
		for (auto&& seq : s)
		{
			std::cout << "> Contig_" << i++ << std::endl;
			std::cout << seq.first /*+ p.second << " " << seq.second << std::endl;
		}
		sequences.insert(sequences.end(), s.begin(), s.end());
	}*/

	std::cerr << (clock() - t)/1000000. << std::endl;
}
