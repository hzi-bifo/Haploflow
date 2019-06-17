#include "Vertex.h"
#include "deBruijnGraph.h"

Vertex::Vertex() : 	
										cc(0),
										index(0),
										onStack(false),
										a_in(0),
										a_out(0),
										c_in(0),
										c_out(0),
										g_in(0),
										g_out(0),
										t_in(0),
										t_out(0),
										terminal(std::make_pair(0,0)),
										n_in(0),
										n_out(0),
										degree(0),
										visited(false),
										flagged(false)
{
}

Vertex::Vertex(int ain, int cin, int gin, int tin, int aout, int cout, int gout, int tout, std::pair<unsigned int, unsigned int> start):
    cc(0),
    index(0),
    onStack(false),
    a_in(ain),
    a_out(aout),
    c_in(cin),
    c_out(cout),
    g_in(gin),
    g_out(gout),
    t_in(tin),
    t_out(tout),
    terminal(start),
    n_in(0),
    n_out(0), // unused right now
    degree(8-!ain-!cin-!gin-!tin-!aout-!cout-!gout-!tout), // degree is 8 if all are set (!x returns 0)
    visited(false),
    flagged(false)
{
}

Vertex::Vertex(const Vertex& v) : 
	cc(v.cc),
	index(v.index),
	onStack(v.onStack),
	a_in(v.a_in),
	a_out(v.a_out),
	c_in(v.c_in),
	c_out(v.c_out),
	g_in(v.g_in),
	g_out(v.g_out),
	t_in(v.t_in),
	t_out(v.t_out),
	terminal(v.terminal),
	n_in(v.n_in),
	n_out(v.n_out),
	degree(v.degree),
	visited(v.visited),
	flagged(v.flagged)
{
}

std::ostream& operator<<(std::ostream& os, const Vertex& v)
{
    os << v.a_in << '\t' << v.c_in << '\t' << v.g_in << '\t' << v.t_in << std::endl;
    os << v.a_out << '\t' << v.c_out << '\t' << v.g_out << '\t' << v.t_out << std::endl;
    os << v.degree << '\t' << v.visited << '\t' << v.terminal.first<< '\t' << v.terminal.second;
    return os;
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
	terminal.first++;
}

void Vertex::read_end()
{
    terminal.second++;
}

void Vertex::visit()
{
	visited = true;
}

void Vertex::unvisit()
{
	visited = false;
}

void Vertex::flag()
{
	flagged = true; // what is flagged may never ... be unflagged
}

bool Vertex::is_flagged() const
{
	return flagged;
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
	return terminal.first;
}

unsigned int Vertex::get_read_ends() const
{
	return terminal.second;
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
