#include <unordered_map>
#include <fstream>
#include <iostream>

class deBruijnGraph
{
public:
	deBruijnGraph(unsigned int k); // creates empty graph
	deBruijnGraph(std::string filename, unsigned int k); // builds the deBruijn graph from file
	void printGraph(); //debug
	
private:
	void add_as_neighbour(const std::string& kmer, const char& letter, bool succ); // add letter as predecessor or successor
	void split_read(const std::string&); // given the read, inserts its kmers in the graph

	std::unordered_map<std::string,std::array<unsigned int,8> > graph_; // successors and predecessors (ACGT)
	unsigned int k_;
};

