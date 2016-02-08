#include <unordered_map>
#include <fstream>

class deBruijnGraph
{
public:
	deBruijnGraph(int k); // creates empty graph
	deBruijnGraph(std::string filename, int k); // builds the deBruijn graph from file
	
private:
	void split_read(std::string); // given the read, inserts its kmers in the graph

	std::unordered_map<std::string,std::array<int,4> > graph_;
	int k_;
};

