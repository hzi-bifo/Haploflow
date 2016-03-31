#include <vector>
#include <iostream>

class Vertex
{
public:
	Vertex();
	Vertex(const Vertex& v);
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
	unsigned int a_in;
	unsigned int a_out;
	unsigned int c_in;
	unsigned int c_out;
	unsigned int g_in;
	unsigned int g_out;
	unsigned int t_in;
	unsigned int t_out;
	unsigned int n_in;
	unsigned int n_out;

};
