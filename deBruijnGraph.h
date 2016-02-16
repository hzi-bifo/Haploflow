#include <unordered_map>
#include <vector>
#include <queue>
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
	std::vector<std::string> get_terminals(bool); // returns all the sinks and/or sources
	void mark_ccs(const std::string&, unsigned int); // calculates pathes from vertex (source) to sinks, returns whether path belongs to a new cc
	void printGraph(); //debug	
	
	std::unordered_map<std::string,Vertex> graph_; // this HAS to be private later on!!!
	// 4 in-edges, 4 out-edges, cc (currently), visited_marker, flow, capacity
	
private:
	void split_read(const std::string&); // given the read, inserts its kmers in the graph
	//std::vector<std::string> get_terminals(bool); // returns all the sinks and/or sources
	
	//std::unordered_map<std::string,std::array<unsigned int,9> > graph_; // successors and predecessors (ACGT) + "flow"
	unsigned int k_;
};

