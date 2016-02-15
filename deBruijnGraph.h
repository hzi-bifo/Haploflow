#include <unordered_map>
#include <vector>
#include <queue>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <limits>

class deBruijnGraph
{
public:
	deBruijnGraph(unsigned int k); // creates empty graph
	deBruijnGraph(std::string filename, unsigned int k); // builds the deBruijn graph from file
	std::vector<std::string> get_terminals(bool); // returns all the sinks and/or sources
	bool mark_ccs(const std::string&, unsigned int, bool); // calculates pathes from vertex (source) to sinks, returns whether path belongs to a new cc
	void printGraph(); //debug
	std::string extractSequence(const std::string&, int); // given a cc (+count of contig)
	
	
	std::unordered_map<std::string,std::array<unsigned int,12> > graph_; // this HAS to be private later on!!!
	// 4 in-edges, 4 out-edges, cc (currently), visited_marker, flow, capacity
	
private:
	void add_as_neighbour(const std::string& kmer, const char& letter, bool succ); // add letter as predecessor or successor
	void split_read(const std::string&); // given the read, inserts its kmers in the graph
	void add_front(int, std::string&);
	void add_back(int, std::string&);
	void edmonds_karp_single(const std::string&, const std::string&);
	// subfunctions of extractSequence depending on the vertex given
	std::string search(const std::string&, bool);
	void resetVisits(); // resets the "visited" bit for another bfs iteration
	//void bfs(const std::string&, int, bool); // calculates pathes from vertex (source) to sinks
	//std::vector<std::string> get_terminals(bool); // returns all the sinks and/or sources
	
	//std::unordered_map<std::string,std::array<unsigned int,9> > graph_; // successors and predecessors (ACGT) + "flow"
	unsigned int k_;
};

