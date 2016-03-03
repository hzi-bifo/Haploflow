#include <vector>
#include <iostream>

class Vertex
{
public:
	Vertex();
	Vertex(const Vertex& v);
	void add_successor(char letter);
	void add_predecessor(char letter);
	const std::vector<std::pair<char, unsigned int> > get_successors() const;
	const std::vector<char> get_predecessors() const;
	const int capacity() const;
	const bool isSource() const;
	const bool isSink() const;
	const bool isJunction() const;
	const void print() const; //debug
	
	unsigned int cc; //these can be set from the outside - handle with care
	bool visited;
	std::string const * source; //points to const string source (of junction)
	
private:
	int a_in;
	int a_out;
	int c_in;
	int c_out;
	int g_in;
	int g_out;
	int t_in;
	int t_out;
	int n_in;
	int n_out;
	int flow;

};
