#ifndef R_INDEX_PHI_SUX_HPP_
#define R_INDEX_PHI_SUX_HPP_

#include <elias_fano_sux.hpp>
#include <sdsl/int_vector.hpp>

namespace stpd{

class r_index_phi_inv_sux
{
public:
	
	r_index_phi_inv_sux(){} // empty constructor

	void build(const std::string bwt_filename, const std::string sa_filename, bool_t verbose = true)
	{
		std::ifstream bwt(bwt_filename,std::ifstream::binary);
		if(not bwt){ std::cerr << "Error opening the BWT file..." << std::endl; exit(1); }
		std::ifstream sa(sa_filename,std::ifstream::binary);
		if(not sa){ std::cerr << "Error opening the SA file..." << std::endl; exit(1); }

		char prev, curr; // previous and current BWT character
		usafe_t curr_sa, prev_sa, bwt_length; // previous and current SA entry
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

		// sort end-of-run samples in ascending order
		std::sort(last_first.begin(), last_first.end(), [](auto &left, auto &right) {
		    return left.first < right.first;
		});

		std::vector<uint64_t> onset; onset.reserve(last_first.size());
		for(auto& idx: last_first){ onset.push_back(idx.first); }
		last.build(onset,bwt_length);

		first = sdsl::int_vector<>(last_first.size(),0,bitsize(uint64_t(bwt_length)));
		if(verbose)
			std::cout << "Number of BWT runs indexed = " << first.size()+1 << std::endl;

		i=0;
		for(auto& entry: last_first)
		{
			first[i++] = entry.second;
		}

		bwt.close();
		sa.close();
	}

	int_t phi_safe(const uint_t idx) const
	{
		if(idx != L)
		{
			auto res = last.successor_rank(idx);
			//std::cout << "-->" << res.first << " - " << res.second << std::endl;
			return first[res.second-1] - (res.first - idx);
		}
		else{ return -1; }
	}

	int_t phi_unsafe(const uint_t idx) const
	{
		auto res = last.successor_rank(idx);

		return first[res.second-1] - (res.first - idx);
	}

	void load(std::istream& in)
	{
		in.read((char*)&L, sizeof(L));
		last.load(in);
		first.load(in);
	}

	uint_t serialize(std::ostream& out)
	{
		uint_t w_bytes = 0;

		out.write((char*)&L, sizeof(L));

		w_bytes += sizeof(L);

		w_bytes += last.serialize(out);
		w_bytes += first.serialize(out);

		return w_bytes;
	}

private:
	
	sux::bits::EliasFano<> last; // last samples
	uint_t L; // last SA entry
	sdsl::int_vector<> first; // first samples
};
}

#endif // R_INDEX_PHI_SUX_HPP_