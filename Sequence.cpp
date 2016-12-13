#include "Sequence.h"

Sequence::Sequence(const std::string& kmer) : 
kmer_(kmer)
{
}

bool Sequence::operator==(const Sequence& seq) const 
{
	return (get_kmer() == seq.get_kmer() or rc() == seq.get_kmer());
}

bool Sequence::operator==(const std::string& s) const
{
	return (get_kmer() == s); // can be used for checking reverse complimentarity
}

bool Sequence::operator!=(const Sequence& seq) const
{
	return !(*this == seq);
}

bool Sequence::operator!=(const std::string& s) const 
{
	return !(*this == s); // BE WARY THAT THIS IS NOT THE SAME EQUALITY LIKE ON THE OBJECT!
}

std::string Sequence::rc() const
{
	std::string rev_compl;
	for (auto&& c : kmer_)
		rev_compl.push_back(complement(c));
	std::reverse(rev_compl.begin(), rev_compl.end());
	return rev_compl;
}

const std::string Sequence::get_kmer() const
{
	return kmer_;
}
