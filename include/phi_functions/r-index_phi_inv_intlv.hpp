#ifndef R_INDEX_PHI_SUX_INTLV_HPP_
#define R_INDEX_PHI_SUX_INTLV_HPP_

#include <elias_fano_intlv.hpp>
#include <common.hpp>

namespace stpd{

class r_index_phi_inv_intlv
{
public:
	
	r_index_phi_inv_intlv(){} // empty constructor

	void build(const std::string bwt_filename, const std::string sa_filename, bool_t verbose = true)
	{
		std::ifstream bwt(bwt_filename,std::ifstream::binary);
		if(not bwt){ std::cerr << "Error opening the BWT file..." << std::endl; exit(1); }
		std::ifstream sa(sa_filename,std::ifstream::binary);
		if(not sa){ std::cerr << "Error opening the SA file..." << std::endl; exit(1); }

		char prev = 0, curr = 0; // previous and current BWT character
		usafe_t curr_sa = 0, prev_sa = 0, bwt_length = 0; // previous and current SA entry
		// Vector of pairs associating each end-of-run sample with
		// the corresponding beginning-of-run sample.
		std::vector<std::pair<usafe_t,usafe_t>> last_first;

		bwt.seekg(0, bwt.end);
		bwt_length = bwt.tellg(); 
		// skip first entry for $
		bwt.seekg(1, bwt.beg); sa.seekg(STORE_SIZE, sa.beg);

		// read the first bwt character and SA entry 
		bwt.read(&prev,1);
		sa.read(reinterpret_cast<char*>(&prev_sa), STORE_SIZE);
		
		usafe_t i=0;
		while(i++<bwt_length)
		{
			bwt.read(&curr,1);
			sa.read(reinterpret_cast<char*>(&curr_sa), STORE_SIZE);

			if(curr != prev)
			{
				last_first.push_back(std::make_pair(prev_sa-1,curr_sa-1));

				prev = curr;
			}

			prev_sa = curr_sa;
		}
		// set last SA entry
		L = prev_sa-1;

		// sort end-of-run samples in increasing order
		std::sort(last_first.begin(), last_first.end(), [](auto &left, auto &right) {
		    return left.first < right.first;
		});

		// construct a sorted dictionary storing (end of run, beginning of run) sample pairs
		LFsamples.build(last_first,bwt_length,bitsize(uint64_t(bwt_length)));

		bwt.close();
		sa.close();
	}

	int_t phi_safe(const uint_t idx) const
	{
		if(idx != L)
		{
			auto res = LFsamples.successor_value(idx);
			return res.second - (res.first - idx);
		}
		else{ return -1; }
	}

	int_t phi_unsafe(const uint_t idx) const
	{
		auto res = LFsamples.successor_value(idx);

		return res.second - (res.first - idx);
	}

	void load(std::istream& in)
	{
		in.read((char*)&L, sizeof(L));

		LFsamples.load(in);
	}

	uint_t serialize(std::ostream& out)
	{
		uint_t w_bytes = 0;

		out.write((char*)&L, sizeof(L));
		w_bytes += sizeof(L);

		w_bytes += LFsamples.serialize(out);

		return w_bytes;
	}

private:

	sux::bits::InterleavedEliasFano<> LFsamples; // last - first samples dictionary
	uint_t L; // last SA entry
};
}

#endif // R_INDEX_PHI_INV_INTLV_HPP_