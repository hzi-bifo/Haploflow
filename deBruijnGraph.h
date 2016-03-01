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
	std::vector<std::string> getScaffolds(std::unordered_map<std::string, std::string>&);
	std::string getSequence(const std::pair<std::string,std::string>&); // calculates the sequence between a source and sink pair
	std::unordered_map<std::string,std::string> find_all_junctions(); // returns pairs of <junction_source,junction_end>
	void printGraph(); //debug	
	int getSize(); // returns number of vertices in graph
	void debug(); // debug
	
	
private:
	std::unordered_map<std::string,Vertex> graph_; //graph data structure
	void split_read(const std::string&); // given the read, inserts its kmers in the graph
	std::string find_next_junction(const std::string&); // find next junction in the graph
	
	unsigned int k_; //kmer size
};

