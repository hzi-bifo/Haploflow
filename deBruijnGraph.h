#ifndef DBG_H
#define DBG_H

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
#include <functional>
#include <cstring>
#include <ctime> //debug
#include "Vertex.h"
#include "Sequence.h"

class deBruijnGraph
{
friend Vertex;
public:
	deBruijnGraph(unsigned int k); // creates empty graph
	deBruijnGraph(std::string filename, unsigned int k); // builds the deBruijn graph from file
	inline static char complement(const char& c){switch (c){case 'A' : return 'T'; case 'C' : return 'G'; case 'G' : return 'C'; case 'T' : return 'A'; default: return 'N';};}
	
	std::vector<std::string> getSources() const; // returns all sources
	std::vector<std::string> getSinks() const; // returns all sinks in graph
	std::pair<std::vector<Sequence>, std::vector<Sequence> > getJunctions() const;
	const Sequence& getSequence(const std::string&);
	Vertex* getVertex(const std::string&); // return the vertex corresponding to kmer, 0 if not found

	void printGraph() const; //debug	
	int getSize() const; // returns number of vertices in graph
	void debug(); // debug
	
private:
	unsigned int split_read(const std::string&); // given the read, inserts its kmers in the graph

	std::unordered_map<Sequence,Vertex> graph_;
	
	unsigned int k_; //kmer size
};

#endif
