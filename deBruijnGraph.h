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
#include "Sequence.h"

class deBruijnGraph
{
public:
	deBruijnGraph(unsigned int k, std::unordered_map<Sequence, Vertex>); // creates graph from subgraph
    deBruijnGraph(std::string filename); // create from file
	deBruijnGraph(std::string filename, unsigned int k); // builds the deBruijn graph from file
	inline static char complement(const char& c){switch (c){case 'A' : return 'T'; case 'C' : return 'G'; case 'G' : return 'C'; case 'T' : return 'A'; default: return 'N';};}
    friend std::ostream& operator<<(std::ostream& os, const deBruijnGraph& dbg);

	void markCycles();
    std::vector<deBruijnGraph> split_ccs();
	std::map<unsigned int, unsigned int> coverageDistribution() const;
	
	std::vector<std::string> getSources() const; // returns all sources
	std::vector<std::string> getSinks() const; // returns all sinks in graph
	std::pair<std::vector<Sequence>, std::vector<Sequence> > getJunctions() const;
	const Sequence* getSequence(const std::string&);
	Vertex* getVertex(const std::string&); // return the vertex corresponding to kmer, 0 if not found

	unsigned int getK() const; // return the value of k
	void printGraph() const; //debug	
	int getSize() const; // returns number of vertices in graph
	void debug(); // debug
	
private:
	unsigned int split_read(const std::string&); // given the read, inserts its kmers in the graph
    std::unordered_map<Sequence, Vertex> dfs(std::pair<Sequence,Vertex>); // depth first search

	std::unordered_map<Sequence, Vertex> graph_;
	
	unsigned int k_; //kmer size
	unsigned int read_length_; // length of longest read
};

#endif
