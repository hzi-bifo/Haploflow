#include <vector>
#include <iostream>
#include <algorithm>

class Vertex
{
	friend bool operator==(const Vertex& lhs, const Vertex& rhs){return (lhs.kmer == rhs.kmer or lhs.rc() == rhs.kmer);};
public:
	Vertex();
	Vertex(const std::string&);
	Vertex(const Vertex& v);
	
	inline bool operator==(const Vertex& rhs){return (kmer == rhs.kmer or rc() == rhs.kmer);};  // these cases are sufficient as rev_compl == rev_compl iff kmer == kmer
	bool isRC(const std::string&) const;
	
	void add_successor(const char& letter) const;
	void add_predecessor(const char& letter) const;
	const std::vector<char> get_successors(bool rc) const;
	const std::vector<char> get_predecessors(bool rc) const; 
	
	const unsigned int capacity(bool rc) const;
	const unsigned int flow(bool rc) const;
	void add_flow(bool rc, unsigned int max_flow) const;
	const bool isSource (bool rc) const;
	const bool isSink(bool rc) const;
	const bool isJunction(bool rc) const;
	const void print(bool cerr) const; //debug
	
	std::string kmer;
	mutable unsigned int cc; //these can be set from the outside - handle with care
	mutable bool visited;
	mutable char pred; // from where the path came
	mutable std::string const* target; //points to const string source (of junction)
	const std::string rc() const;
	
	
private:
	inline static char complement(const char& c){switch (c){case 'A' : return 'T'; case 'C' : return 'G'; case 'G' : return 'C'; case 'T' : return 'A'; default: return 'N';};}
	//std::string rev_compl;
	mutable unsigned int flow_f, flow_r;
	mutable unsigned int a_in, a_in_r;
	mutable unsigned int a_out, a_out_r;
	mutable unsigned int c_in, c_in_r;
	mutable unsigned int c_out, c_out_r;
	mutable unsigned int g_in, g_in_r;
	mutable unsigned int g_out, g_out_r;
	mutable unsigned int t_in, t_in_r;
	mutable unsigned int t_out, t_out_r;
	mutable unsigned int n_in, n_in_r;
	mutable unsigned int n_out, n_out_r;

};

namespace std {
  template <> struct hash<Vertex>
  {
    size_t operator()(const Vertex& x) const
    {
    	auto v = hash<string>()(x.kmer);
    	auto w = hash<string>()(x.rc());
      return min(v,w);
    }
  };
}
