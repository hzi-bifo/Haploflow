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

class deBruijnGraph
{
public:
	deBruijnGraph(unsigned int k); // creates empty graph
	deBruijnGraph(std::string filename, bool fasta, unsigned int k); // builds the deBruijn graph from file
	
	std::pair<std::string,unsigned int> bfs(const std::string&, int*, std::function<void(const std::string&, int*)>, std::function<bool(const std::string&, int*)>, bool); // generalized bfs
	void add_sequence(std::string filename);

	std::vector<std::string> getSources(); // returns all sources
	std::vector<std::string> getSinks(); // returns all sinks in graph
	std::unordered_map<std::string,std::string> find_all_junctions(); // returns pairs of <junction_source,junction_end>
	std::vector<std::pair<std::string, unsigned int> > getSequences(const std::string& source, const std::string& sink); // finds all paths + sequences between two junctions
	std::pair<std::string,unsigned int> glue(const std::string& source, const std::unordered_map<std::string, std::string>&); // glues together possible contigs
	
	std::string make_graph();	

	void printGraph(); //debug	
	int getSize(); // returns number of vertices in graph
	void debug(); // debug
	
	
private:
	inline static char complement(const char& c){switch (c){case 'A' : return 'T'; case 'C' : return 'G'; case 'G' : return 'C'; case 'T' : return 'A'; default: return 'N';};}
	unsigned int split_read(const std::string&); // given the read, inserts its kmers in the graph
	std::string find_next_junction(const std::string*); // find next junction in the graph
	void split_fasta(std::string filename); // split fasta into kmers

	std::unordered_set<Vertex> graph_;
	unsigned int k_; //kmer size
};
