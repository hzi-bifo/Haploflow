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
	void split_read(std::string); // given the read, inserts its kmers in the graph

	std::unordered_map<std::string,std::array<unsigned int,4> > graph_; // ACGT
	unsigned int k_;
};

