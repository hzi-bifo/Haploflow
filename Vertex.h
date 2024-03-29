#ifndef VERTEX_H
#define VERTEX_H

#include <vector>
#include <iostream>

class deBruijnGraph;
class Sequence; 

class Vertex
{
public:
    Vertex();
	Vertex(int ain, int cin, int gin, int tin, int aout, int cout, int gout, int tout, std::pair<unsigned int, unsigned int> terminal);
	Vertex(const Vertex& v);
    friend std::ostream& operator<<(std::ostream& os, const Vertex& v);
	
	void add_successor(const char& letter);
	void add_predecessor(const char& letter);
	void read_start();
	void read_end();
	void visit();
	void unvisit();
	void flag();
	
	unsigned int cc; // connected component
	unsigned int index; // to possibly identify this vertex in the unitig graph
	bool onStack; // dfs-specific variable

	bool is_junction() const;
	bool is_conflicting() const;
	bool is_visited() const;
	bool is_flagged() const;
	std::vector<char> get_successors() const;
	std::vector<char> get_predecessors() const; 
	unsigned int get_read_starts() const;
	unsigned int get_read_ends() const;
	unsigned int get_out_coverage(const char) const;
	unsigned int get_in_coverage(const char) const;
	unsigned int get_total_in_coverage() const;
	unsigned int get_total_out_coverage() const;
	unsigned int get_degree() const;

	const bool isSource () const;
	const bool isSink() const;
	const void print(bool cerr) const; //debug

private:
	unsigned int a_in;
	unsigned int a_out;
	unsigned int c_in;
	unsigned int c_out;
	unsigned int g_in;
	unsigned int g_out;
	unsigned int t_in;
	unsigned int t_out;
	std::pair<unsigned int, unsigned int> terminal; // number of reads starting with *this
	unsigned int n_in;
	unsigned int n_out;
	unsigned int degree; // in_degree + out_degree
	bool visited;
	bool flagged; // flagged for delete
};

#endif
