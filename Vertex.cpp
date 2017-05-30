#include "Vertex.h"
#include "deBruijnGraph.h"

Vertex::Vertex() : 	
										a_in(0),
										a_out(0),
										c_in(0),
										c_out(0),
										g_in(0),
										g_out(0),
										t_in(0),
										t_out(0),
										n_in(0),
										n_out(0),
										degree(0),
										visited(false),
										flagged(false),
										index(0),
										starts_with(0),
										cycle(false)
{
}

Vertex::Vertex(const Vertex& v) : 
	a_in(v.a_in),
	a_out(v.a_out),
	c_in(v.c_in),
	c_out(v.c_out),
	g_in(v.g_in),
	g_out(v.g_out),
	t_in(v.t_in),
	t_out(v.t_out),
	n_in(v.n_in),
	n_out(v.n_out),
	degree(v.degree),
	visited(v.visited),
	flagged(v.flagged),
	index(v.index),
	starts_with(v.starts_with),
	cycle(v.cycle)
{
}

// if A -> B with letter C, then rc(B) -> rc(A) with c(C)
void Vertex::add_successor(const char& letter)
{
	bool is_set = false; //check whether this successor already has been set
	switch(letter)
	{
		case 'A': is_set = a_out; a_out++; break;
		case 'C': is_set = c_out; c_out++; break;
		case 'G': is_set = g_out; g_out++; break;
		case 'T': is_set = t_out; t_out++; break;
		default: is_set = n_out; n_out++; break;
	}
	if (!is_set) degree++;
	// junction if it has been visited (has successors) but this particular successor has not been set
	// this means it has at least 2 different successors
}

// the same for predecessors
void Vertex::add_predecessor(const char& letter)
{
	bool is_set = false;
	switch(letter)
	{
		case 'A': is_set = a_in; a_in++; break;
		case 'C': is_set = c_in; c_in++; break;
		case 'G': is_set = g_in; g_in++; break;
		case 'T': is_set = t_in; t_in++; break;
		default: is_set = n_in; n_in++; break;
	}
	if (!is_set) degree++;
}

void Vertex::read_start()
{
	starts_with++;
}

void Vertex::visit()
{
	visited = true;
}

void Vertex::unvisit()
{
	visited = false;
}

void Vertex::mark_as_cycle()
{
	cycle = true;
}

bool Vertex::is_cycle() const
{
	return cycle;
}

void Vertex::flag()
{
	flagged = true; // what is flagged may never ... be unflagged
}

bool Vertex::is_flagged() const
{
	return flagged;
}

void Vertex::set_index(unsigned int i)
{
	index = i;
}

unsigned int Vertex::get_index() const
{
	return index;
}

bool Vertex::is_visited() const
{
	return visited;
}

// if indegree not equal outdegree
bool Vertex::is_junction() const
{
	return (!!a_in + !!c_in + !!g_in + !!t_in + !!n_in != (!!a_out + !!c_out + !!g_out + !!t_out + !!n_out));
}

bool Vertex::is_conflicting() const
{
	return ((!!a_in + !!c_in + !!g_in + !!t_in + !!n_in) > 1
	and (!!a_out + !!c_out + !!g_out + !!t_out + !!n_out) > 1 and !is_junction());
}

unsigned int Vertex::get_degree() const
{
	return degree;
}

// TODO beauty
// returns the forward successors
std::vector<char> Vertex::get_successors() const
{
	std::vector<char> ret;
	if (a_out) ret.push_back('A');
	if (c_out) ret.push_back('C');
	if (g_out) ret.push_back('G');
	if (t_out) ret.push_back('T');
	if (n_out) ret.push_back('N');
	return ret;
}

// returns the forward predecessors
std::vector<char> Vertex::get_predecessors() const
{
	std::vector<char> ret;
	if (a_in) ret.push_back('A');
	if (c_in) ret.push_back('C');
	if (g_in) ret.push_back('G');
	if (t_in) ret.push_back('T');
	if (n_in) ret.push_back('N');
	return ret;
}

unsigned int Vertex::get_read_starts() const
{
	return starts_with;
}

unsigned int Vertex::get_out_coverage(const char c) const
{
	if (c == 'A') return a_out;
	else if (c == 'C') return c_out;
	else if (c == 'G') return g_out;
	else if (c == 'T') return t_out;
	else return n_out;
}

unsigned int Vertex::get_in_coverage(const char c) const
{
	if (c == 'A') return a_in;
	else if (c == 'C') return c_in;
	else if (c == 'G') return g_in;
	else if (c == 'T') return t_in;
	else return n_in;
}

unsigned int Vertex::get_total_out_coverage() const
{
	return a_out + c_out + g_out + t_out + n_out;
}

unsigned int Vertex::get_total_in_coverage() const
{
	return a_in + c_in + g_in + t_in + n_in;
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
