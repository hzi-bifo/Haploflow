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
	
	bool is_junction() const;
	bool is_conflicting() const;
	const std::vector<unsigned int> get_successors() const;
	const std::vector<unsigned int> get_predecessors() const; 
	
	const bool isSource () const;
	const bool isSink() const;
	const void print(bool cerr) const; //debug

private:
	mutable unsigned int a_in;
	mutable unsigned int a_out;
	mutable unsigned int c_in;
	mutable unsigned int c_out;
	mutable unsigned int g_in;
	mutable unsigned int g_out;
	mutable unsigned int t_in;
	mutable unsigned int t_out;
	mutable unsigned int n_in;
	mutable unsigned int n_out;
	bool visited;
};

#endif
