#include "Vertex.h"

Vertex::Vertex() : 	kmer(""),
										//rev_compl(""),
										flow(0),
										visited(false),
										a_in(0), a_in_r(0),
										a_out(0), a_out_r(0),
										c_in(0), c_in_r(0),
										c_out(0), c_out_r(0),
										g_in(0), g_in_r(0),
										g_out(0), g_out_r(0),
										t_in(0), t_in_r(0),
										t_out(0), t_out_r(0),
										n_in(0), n_in_r(0),
										n_out(0), n_out_r(0),
{
}

Vertex::Vertex(const std::string& kmer) : 	
										kmer(kmer),
										//rev_compl(rc(kmer)),
										flow(0),
										visited(false),
										a_in(0), a_in_r(0),
										a_out(0), a_out_r(0),
										c_in(0), c_in_r(0),
										c_out(0), c_out_r(0),
										g_in(0), g_in_r(0),
										g_out(0), g_out_r(0),
										t_in(0), t_in_r(0),
										t_out(0), t_out_r(0),
										n_in(0), n_in_r(0),
										n_out(0), n_out_r(0)
{
}

Vertex::Vertex(const Vertex& v) : 
	kmer(v.kmer),
	//rev_compl(v.rev_compl),
	cc(v.cc),
	flow(v.flow),
	visited(v.visited),
	pred(v.pred),
	source(v.source),
	a_in(v.a_in), a_in_r(v.a_in_r),
	a_out(v.a_out), a_out_r(v.a_out_r),
	c_in(v.c_in), c_in_r(v.c_in_r),
	c_out(v.c_out), c_out_r(v.c_out_r),
	g_in(v.g_in), g_in_r(v.g_in_r),
	g_out(v.g_out), g_out_r(v.g_out_r),
	t_in(v.t_in), t_in_r(v.t_in_r),
	t_out(v.t_out), t_out_r(v.t_out_r),
	n_in(v.n_in), n_in_r(v.n_in_r),
	n_out(v.n_out), n_out_r(v.n_out_r),
{
}

// returns reverse complement of a string
std::string Vertex::rc (const std::string& kmer)
{
	std::string rev(kmer);
	std::transform(kmer.begin(),kmer.end(),rev.begin(),complement);
	std::reverse(rev.begin(),rev.end());
	return rev;
}

bool Vertex::isRC(const std::string& s)
{
	return (rc(s) == kmer);
}

// if A -> B with letter C, then rc(B) -> rc(A) with c(C)
void Vertex::add_successor(char letter)
{
	switch(letter)
	{
		case 'A': a_out++; t_in_r++; break;
		case 'C': c_out++; g_in_r++; break;
		case 'G': g_out++; c_in_r++; break;
		case 'T': t_out++; a_in_r++; break;
		default: n_out++; n_in_r++; break;
	}
	
}

// the same for predecessors
void Vertex::add_predecessor(char letter)
{
	switch(letter)
	{
		case 'A': a_in++; t_out_r++; break;
		case 'C': c_in++; g_out_r++; break;
		case 'G': g_in++; c_out_r++; break;
		case 'T': t_in++; a_out_r++; break;
		default: n_in++; n_out_r++; break;
	}
}

// returns the forward successors
const std::vector<char> Vertex::get_successors() const
{
	std::vector<char> succ;
	if (a_out)
		succ.push_back('A');
	if (c_out)
		succ.push_back('C');
	if (g_out)
		succ.push_back('G');
	if (t_out)
		succ.push_back('T');
	if (n_out)
		succ.push_back('N');
	return succ;
}

// returns the forward predecessors
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

const unsigned int Vertex::capacity() const
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
	std::cout << kmer << "/" << rev_compl << std::endl;
	std::cout << "Out - A: " << a_out << ", C: " << c_out << ", G: " << g_out << ", T: " << t_out << ", N: " << n_out << std::endl;
	std::cout << "In  - A: " << a_in << ", C: " << c_in << ", G: " << g_in << ", T: " << t_in << ", N: " << n_in << std::endl;
	std::cout << "Component: " << cc << ", capacity: " << capacity() << ", used flow: " << flow << std::endl;
}
