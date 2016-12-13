#ifndef DBG_H
#define DBG_H

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
#include "UnitigGraph.h"

class deBruijnGraph
{
public:
	deBruijnGraph(unsigned int k); // creates empty graph
	deBruijnGraph(std::string filename, bool fasta, unsigned int k); // builds the deBruijn graph from file
	
	Vertex find(const std::string&);

	template<typename T>
	std::pair<std::string,unsigned int> bfs(const std::string&, T*, std::function<void(const std::string&, T*)>, std::function<bool(const std::string&, T*)>, bool); // generalized bfs
	template<typename T>
	std::pair<std::string,unsigned int> dfs(const std::string&, T*, std::function<void(const std::string&, T*)>, std::function<bool(const std::string&, T*)>, bool); // generalized bfs
	void add_sequence(std::string filename);

	std::vector<std::string> getSources(); // returns all sources
	std::vector<std::string> getSinks(); // returns all sinks in graph
	std::vector<std::pair<std::string, unsigned int> > getSequences(const std::string& source, const std::string& sink); // finds all paths + sequences between two junctions
	
	std::string make_graph();	

	void printGraph(); //debug	
	int getSize(); // returns number of vertices in graph
	void debug(); // debug
	
private:
	inline static char complement(const char& c){switch (c){case 'A' : return 'T'; case 'C' : return 'G'; case 'G' : return 'C'; case 'T' : return 'A'; default: return 'N';};}
	unsigned int split_read(const std::string&); // given the read, inserts its kmers in the graph
	void split_fasta(std::string filename); // split fasta into kmers
	std::pair<std::vector<Vertex>,std::vector<Vertex> > getUnbalanced(); // returns all unbalanced vertices

	std::unordered_set<Vertex> graph_;
	
	unsigned int k_; //kmer size
};

#endif
