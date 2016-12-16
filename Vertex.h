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
	Vertex(const Vertex& v);
	
	void add_successor(const char& letter);
	void add_predecessor(const char& letter);
	void visit();
	void unvisit();
	void set_index(unsigned int);
	unsigned int get_index() const;

	bool is_junction() const;
	bool is_conflicting() const;
	bool is_visited() const;
	std::vector<char> get_successors() const;
	std::vector<char> get_predecessors() const; 
	
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
	unsigned int n_in;
	unsigned int n_out;
	bool visited;
	unsigned int index; // to possibly identify this vertex in the unitig graph
};

#endif
