#include <unordered_map>
#include <unordered_set>
#include <map>
#include <vector>
#include <queue>
#include <stack>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <functional>
#include <cstring>
#include <ctime> //debug
#include "Vertex.h"

class deBruijnGraph
{
public:
	deBruijnGraph(unsigned int k); // creates empty graph
	deBruijnGraph(std::string filename, unsigned int k); // builds the deBruijn graph from file
	
	std::pair<std::string,unsigned int> bfs(const std::string&, bool, std::function<void(std::string&)>, std::function<bool(std::string&)>, bool); // generalized bfs
	std::string dfs(const std::string&, bool, std::function<void(std::string&)>, std::function<bool(std::string&)>, bool); // generalized dfs

	std::vector<std::string> getSources(); // returns all sources
	std::vector<std::string> getSinks(); // returns all sinks in graph
	std::unordered_map<std::string,std::string> find_all_junctions(); // returns pairs of <junction_source,junction_end>
	std::vector<std::pair<std::string, unsigned int> > getSequences(const std::string& source, const std::string& sink); // finds all paths + sequences between two junctions
	void printGraph(); //debug	
	int getSize(); // returns number of vertices in graph
	void debug(); // debug
	
	
private:
	void split_read(const std::string&); // given the read, inserts its kmers in the graph
	std::string find_next_junction(const std::string*); // find next junction in the graph
	std::string reverse_complement(const std::string&);
	
	std::function<char(const char&)> lambda = [](const char& c){switch (c){case 'A' : return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T' : return 'A'; default : return 'N';}};
	
	std::unordered_map<std::string,Vertex/*,SeqHash,SeqEq*/> graph_; //graph data structure
	unsigned int k_; //kmer size
};


struct SeqHash
{
  unsigned long long operator() (const std::string& s) const
  {
  	size_t i = 0;
  	unsigned long long hash = 0;
  	unsigned long long rchash = 0;
  	auto lambda = [](const char& c){ switch (c){case 'A' : return 0; case 'C' : return 1; case 'G': return 2; case 'T': return 3; default: return 4;}; }; //0b4?
  	const char* d1 = s.data();
  	while (std::strlen(d1))
  	{
  		rchash += (3 - lambda(*(d1))) * std::pow(4,i++);
  		hash *= 4; 
  		hash += lambda(*(d1++));
  	}
  	return std::min(hash,rchash);
  }
};

struct SeqEq
{
	bool operator() (const std::string& s1, const std::string& s2) const
	{
		std::string rev(s1);
		std::transform(s1.begin(),s1.end(),rev.begin(),[](const char& c){switch (c){case 'A' : return 'T'; case 'C' : return 'G'; case 'G' : return 'C'; case 'T' : return 'A'; default: return 'N';};});
		std::reverse(rev.begin(),rev.end());
		return (s1 == s2 or rev == s2);
	}
};
