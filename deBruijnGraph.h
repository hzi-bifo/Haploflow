#include <unordered_map>
#include <vector>
#include <queue>
#include <fstream>
#include <iostream>

class deBruijnGraph
{
public:
	deBruijnGraph(unsigned int k); // creates empty graph
	deBruijnGraph(std::string filename, unsigned int k); // builds the deBruijn graph from file
	std::vector<std::string> get_terminals(bool); // returns all the sinks and/or sources
	void printGraph(); //debug
	
private:
	void add_as_neighbour(const std::string& kmer, const char& letter, bool succ); // add letter as predecessor or successor
	void split_read(const std::string&); // given the read, inserts its kmers in the graph
	void bfs(const std::string&, int state); // calculates pathes from vertex (source) to sinks + assigns flow
	void backtracking(const std::string&);
	
	std::unordered_map<std::string,std::array<unsigned int,9> > graph_; // successors and predecessors (ACGT) + "flow"
	unsigned int k_;
};

