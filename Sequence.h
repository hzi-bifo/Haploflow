#ifndef SEQ
#define SEQ

#include <string>
#include <algorithm>
#include "nthash.hpp"

class Sequence
{
public:
	Sequence(const std::string& kmer);
	bool operator==(const Sequence& seq) const;
	bool operator==(const std::string& s) const;
	bool operator!=(const Sequence& seq) const;
	bool operator!=(const std::string& s) const;
	friend std::ostream& operator<<(std::ostream& os, const Sequence& s){return os << s.get_kmer();};
	const std::string get_kmer() const;
    std::string rc() const;
private:
	inline static char complement(char c){switch(c){ case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T' : return 'A'; default : return 'N';};}
	std::string kmer_;
};


namespace std 
{
	template <>
	struct hash<Sequence>
	{
		size_t operator()(const Sequence& k) const
		{
            std::string toHash = k.get_kmer();
            return NTC64(toHash.c_str(), toHash.size());
		}
	};
}

#endif
