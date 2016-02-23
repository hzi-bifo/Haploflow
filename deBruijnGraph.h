#include <unordered_map>
#include <map>
#include <vector>
#include <queue>
#include <stack>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <limits>
#include "Vertex.h"

class deBruijnGraph
{
public:
	deBruijnGraph(unsigned int k); // creates empty graph
	deBruijnGraph(std::string filename, unsigned int k); // builds the deBruijn graph from file
	std::string bfs(const std::string&, int, bool, void (*)(Vertex&, int), bool (*)(const Vertex&, int), bool); // generalized bfs
	std::string dfs(const std::string&, int, bool, void (*)(Vertex&, int), bool (*)(const Vertex&, int), bool); // generalized dfs
	std::vector<std::string> getSources(); // returns all sources
	std::vector<std::string> getSinks(); // returns all sinks in graph
	void printGraph(); //debug	
	int getSize(); // returns number of vertices in graph
	
//going to be private:	
	std::string find_next_junction(const std::string&); // find next junction in the graph
	
private:
	std::unordered_map<std::string,Vertex> graph_; //graph data structure
	void split_read(const std::string&); // given the read, inserts its kmers in the graph
	//std::string find_next_junction(const std::string&); // find next junction in the graph
	
	//std::unordered_map<std::string,std::array<unsigned int,9> > graph_; // successors and predecessors (ACGT) + "flow"
	unsigned int k_;
};

