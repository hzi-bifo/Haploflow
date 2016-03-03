#include "Vertex.h"

Vertex::Vertex() : 	cc(0),
										visited(false),
										source(0),
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
										flow(0)
{
}

Vertex::Vertex(const Vertex& v) : 
	cc(v.cc),
	visited(v.visited),
	source(v.source),
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
	flow(v.flow)
{
}

void Vertex::add_successor(char letter)
{
	switch(letter)
	{
		case 'A': a_out++; break;
		case 'C': c_out++; break;
		case 'G': g_out++; break;
		case 'T': t_out++; break;
		default: n_out++; break;
	}
}

void Vertex::add_predecessor(char letter)
{
	switch(letter)
	{
		case 'A': a_in++; break;
		case 'C': c_in++; break;
		case 'G': g_in++; break;
		case 'T': t_in++; break;
		default: n_in++; break;
	}
}

const std::vector<std::pair<char, unsigned int> > Vertex::get_successors() const
{
	std::vector<std::pair<char, unsigned int> > succ;
	if (a_out)
		succ.push_back(std::make_pair('A',a_out));
	if (c_out)
		succ.push_back(std::make_pair('C',c_out));
	if (g_out)
		succ.push_back(std::make_pair('G',g_out));
	if (t_out)
		succ.push_back(std::make_pair('T',t_out));
	if (n_out)
		succ.push_back(std::make_pair('N',n_out));
	return succ;
}

const std::vector<char> Vertex::get_predecessors() const
{
	std::vector<char> pred;
	if (a_in)
		pred.push_back('A');
	if (c_in)
		pred.push_back('C');
	if (g_in)
		pred.push_back('G');
	if (t_in)
		pred.push_back('T');
	if (n_in)
		pred.push_back('N');
	return pred;
}

const int Vertex::capacity() const
{
	int cap1 = a_out + c_out + g_out + t_out;
	return cap1; // this is the "out-capacity"
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

const bool Vertex::isJunction() const
{
	bool diff1 = ((((a_out xor c_out) xor g_out) xor t_out) xor n_out);
	bool diff2 = ((((a_in xor c_in) xor g_in) xor t_in) xor n_in);
	return (diff1 or diff2);
}

//debug
const void Vertex::print() const
{
	std::cout << "Out - A: " << a_out << ", C: " << c_out << ", G: " << g_out << ", T: " << t_out << ", N: " << n_out << std::endl;
	std::cout << "In  - A: " << a_in << ", C: " << c_in << ", G: " << g_in << ", T: " << t_in << ", N: " << n_in << std::endl;
	std::cout << "Component: " << cc << ", capacity: " << capacity() << ", used flow: " << flow << std::endl;
}
