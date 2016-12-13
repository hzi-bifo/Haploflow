#include "Vertex.h"
#include "deBruijnGraph.h"

Vertex::Vertex() : 	
										visited(false),
										a_in(0),
										a_out(0),
										c_in(0),
										c_out(0),
										g_in(0),
										g_out(0),
										t_in(0),
										t_out(0),
										n_in(0),
										n_out(0)
{
}

Vertex::Vertex(const Vertex& v) : 
	visited(v.visited),
	a_in(v.a_in),
	a_out(v.a_out),
	c_in(v.c_in),
	c_out(v.c_out),
	g_in(v.g_in),
	g_out(v.g_out),
	t_in(v.t_in),
	t_out(v.t_out),
	n_in(v.n_in),
	n_out(v.n_out)
{
}

// if A -> B with letter C, then rc(B) -> rc(A) with c(C)
void Vertex::add_successor(const char& letter, deBruijnGraph& dbg)
{
	bool is_zero = true;
	switch(letter)
	{
		case 'A': is_zero = a_out; a_out++; break;
		case 'C': is_zero = c_out; c_out++; break;
		case 'G': is_zero = g_out; g_out++; break;
		case 'T': is_zero = t_out; t_out++; break;
		default: is_zero = n_out; n_out++; break;
	}
}

void Vertex::add_successor(Vertex* next)
{
	successors.push_back(next);
}

// the same for predecessors
void Vertex::add_predecessor(const char& letter, deBruijnGraph& dbg)
{
	bool is_zero = true;
	switch(letter)
	{
		case 'A': is_zero = a_in; a_in++; break;
		case 'C': is_zero = c_in; c_in++; break;
		case 'G': is_zero = g_in; g_in++; break;
		case 'T': is_zero = t_in; t_in++; break;
		default: is_zero = n_in; n_in++; break;
	}
}

void Vertex::add_predecessor(Vertex* prev)
{
	predecessors.push_back(prev);
}

// returns the forward successors
const std::vector<Vertex*> Vertex::get_successors() const
{
	return successors;
}

// returns the forward predecessors
const std::vector<Vertex*> Vertex::get_predecessors() const
{
	return predecessors;
}

const bool Vertex::isSource() const
{
	bool source = !(a_in + c_in + g_in + t_in + n_in); // all 0 ("real" source)
	return source;
}

const bool Vertex::isSink() const
{
	bool sink = !(a_out + c_out + g_out + t_out + n_out); // all 0 ("real" sink)
	return sink;
}

//debug
const void Vertex::print(bool cerr) const
{
	if (cerr)
	{
		std::cerr << "Out - A: " << a_out << ", C: " << c_out << ", G: " << g_out << ", T: " << t_out << ", N: " << n_out << std::endl;
		std::cerr << "In  - A: " << a_in << ", C: " << c_in << ", G: " << g_in << ", T: " << t_in << ", N: " << n_in << std::endl;
	}
	else
	{
		std::cout << "Out - A: " << a_out << ", C: " << c_out << ", G: " << g_out << ", T: " << t_out << ", N: " << n_out << std::endl;
		std::cout << "In  - A: " << a_in << ", C: " << c_in << ", G: " << g_in << ", T: " << t_in << ", N: " << n_in << std::endl;
	}
}
