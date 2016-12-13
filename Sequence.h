#include <algorithm>

class Sequence
{
public:
	Sequence(const std::string& kmer);
	bool operator==(const Sequence& seq) const;
	bool operator==(const std::string& s) const;
	bool operator!=(const Sequence& seq) const;
	bool operator!=(const std::string& s) const;
	std::string rc() const;
	const std::string get_kmer() const;
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
			return (min(hash<string>()(k.get_kmer()),hash<string>()(k.rc())));
		}
	};
}
