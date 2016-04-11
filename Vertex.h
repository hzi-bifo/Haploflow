#include <vector>
#include <iostream>

class Vertex
{
public:
	Vertex();
	Vertex(const std::string&);
	Vertex(const Vertex& v);
	
	inline bool operator==(const Vertex& lhs, const Vertex& rhs){return (lhs.kmer == rhs.kmer or lhs.rev_compl == rhs.kmer);} // these cases are sufficient as rev_compl == rev_compl iff kmer == kmer
	bool isRC(const std::string&);
	
	void add_successor(char letter);
	void add_predecessor(char letter);
	const std::vector<char> get_successors() const;
	const std::vector<char> get_predecessors() const;
	const unsigned int capacity() const;
	const bool isSource() const;
	const bool isSink() const;
	const bool isJunction() const;
	const void print() const; //debug
	
	unsigned int cc; //these can be set from the outside - handle with care
	unsigned int flow;
	bool visited;
	char pred; // from where the path came
	std::string const* source; //points to const string source (of junction)
	
private:
	auto complement[](const char& c){switch (c){case 'A' : return 'T'; case 'C' : return 'G'; case 'G' : return 'C'; case 'T' : return 'A'; default: return 'N';};}
	std::string rc (const std::string&);
	std::string kmer;
	//std::string rev_compl;
	unsigned int a_in, a_in_r;
	unsigned int a_out, a_out_r;
	unsigned int c_in, c_in_r;
	unsigned int c_out, c_out_r;
	unsigned int g_in, g_in_r;
	unsigned int g_out, g_out_r;
	unsigned int t_in, t_in_r;
	unsigned int t_out, t_out_r;
	unsigned int n_in, n_in_r;
	unsigned int n_out, n_out_r;

};
