#include "Vertex.h"

Vertex::Vertex() : 	kmer(""),
										cc(0),
										visited(false),
										flow_f(0), flow_r(0),
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

Vertex::Vertex(const std::string& kmer) : 	
										kmer(kmer),
										cc(0),
										visited(false),
										flow_f(0), flow_r(0),
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
	cc(v.cc),
	visited(v.visited),
	pred(v.pred),
	target(v.target),
	flow_f(v.flow_f), flow_r(v.flow_r),
	a_in(v.a_in), a_in_r(v.a_in_r),
	a_out(v.a_out), a_out_r(v.a_out_r),
	c_in(v.c_in), c_in_r(v.c_in_r),
	c_out(v.c_out), c_out_r(v.c_out_r),
	g_in(v.g_in), g_in_r(v.g_in_r),
	g_out(v.g_out), g_out_r(v.g_out_r),
	t_in(v.t_in), t_in_r(v.t_in_r),
	t_out(v.t_out), t_out_r(v.t_out_r),
	n_in(v.n_in), n_in_r(v.n_in_r),
	n_out(v.n_out), n_out_r(v.n_out_r)
{
}

// returns reverse complement of a string
const std::string Vertex::rc() const
{
	std::string rev(kmer);
	std::transform(kmer.begin(),kmer.end(),rev.begin(),complement);
	std::reverse(rev.begin(),rev.end());
	return rev;
}

bool Vertex::isRC(const std::string& s) const
{
	return (s != kmer);
}

// if A -> B with letter C, then rc(B) -> rc(A) with c(C)
void Vertex::add_successor(const char& letter) const
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
void Vertex::add_predecessor(const char& letter) const
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
const std::vector<char> Vertex::get_successors(bool rc) const
{
	std::vector<char> succ;
	if ((!rc and a_out) or (a_out_r and rc))
		succ.push_back('A');
	if ((!rc and c_out) or (c_out_r and rc))
		succ.push_back('C');
	if ((!rc and g_out) or (g_out_r and rc))
		succ.push_back('G');
	if ((!rc and t_out) or (t_out_r and rc))
		succ.push_back('T');
	if ((!rc and n_out) or (n_out_r and rc))
		succ.push_back('N');
	return succ;
}

// returns the forward predecessors
const std::vector<char> Vertex::get_predecessors(bool rc) const
{
	std::vector<char> pred;
	if ((!rc and a_in) or (a_in_r and rc))
		pred.push_back('A');
	if ((!rc and c_in) or (c_in_r and rc))
		pred.push_back('C');
	if ((!rc and g_in) or (g_in_r and rc))
		pred.push_back('G');
	if ((!rc and t_in) or (t_in_r and rc))
		pred.push_back('T');
	if ((!rc and n_in) or (n_in_r and rc))
		pred.push_back('N');
	return pred;
}

const unsigned int Vertex::capacity(bool rc) const
{
	int cap1 = a_in + c_in + g_in + t_in;
	int cap2 = a_in_r + c_in_r + g_in_r + t_in_r;
	return (rc ? cap2 : cap1); // this is the "in-capacity". We do not check capacity of source, only of the following!
}

const unsigned int Vertex::flow(bool rc) const
{
	return (rc ? flow_r : flow_f);
}

void Vertex::add_flow(bool rc, unsigned int max_flow) const
{
	(rc ? flow_r : flow_f) += max_flow;
}
const bool Vertex::isSource(bool rc) const
{
	bool source = !(a_in + c_in + g_in + t_in + n_in); // all 0 ("real" source)
	return ((!rc and source) or (rc and isSink(!rc)));
}

const bool Vertex::isSink(bool rc) const
{
	bool sink = !(a_out + c_out + g_out + t_out + n_out); // all 0 ("real" sink)
	return ((!rc and sink) or (rc and isSource(!rc))); // if we are reverse complement and sink, then original is a source 
}

const bool Vertex::isJunction(bool rc) const
{
	bool diff1 = ((((a_out xor c_out) xor g_out) xor t_out) xor n_out);
	bool diff2 = ((((a_in xor c_in) xor g_in) xor t_in) xor n_in);
	return (diff1 or diff2);
}

//debug
const void Vertex::print(bool cerr) const
{
	if (cerr)
	{
		std::cerr << kmer << "/" << rc() << std::endl;
		std::cerr << "Out - A: " << a_out << ", C: " << c_out << ", G: " << g_out << ", T: " << t_out << ", N: " << n_out << std::endl;
		std::cerr << "In  - A: " << a_in << ", C: " << c_in << ", G: " << g_in << ", T: " << t_in << ", N: " << n_in << std::endl;
		std::cerr << "Component: " << cc << ", capacity: " << capacity(false) << ", used flow: " << flow_f << " (forward)/ " << flow_r << " (reverse)" << std::endl;
	}
	else
	{
		std::cout << kmer << "/" << rc() << std::endl;
		std::cout << "Out - A: " << a_out << ", C: " << c_out << ", G: " << g_out << ", T: " << t_out << ", N: " << n_out << std::endl;
		std::cout << "In  - A: " << a_in << ", C: " << c_in << ", G: " << g_in << ", T: " << t_in << ", N: " << n_in << std::endl;
		std::cout << "Component: " << cc << ", capacity: " << capacity(false) << ", used flow: " << flow_f << " (forward)/ " << flow_r << " (reverse)" << std::endl;
	}
}
