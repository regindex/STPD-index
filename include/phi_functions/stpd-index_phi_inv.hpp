#ifndef STPD_INDEX_PHI_HPP_
#define STPD_INDEX_PHI_HPP_

#include <elias_fano_bitvector.hpp>
#include <sdsl/int_vector.hpp>

namespace stpd{

template<class b_v = bv::elias_fano_bitvector>
class stpd_index_phi_inv
{
public:
	// empty constructor
	stpd_index_phi_inv(){}

	void build(const std::string bwt_filename, const std::string sa_filename,
		       const sdsl::int_vector<>* stpdA_, usafe_t N)
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

		bwt.close();
		sa.close();

		std::cout << "LAst first vector size= " << last_first.size() << std::endl;

		// sort end-of-run samples in ascending order
		std::sort(last_first.begin(), last_first.end(), [](auto &left, auto &right) {
		    return left.first < right.first;
		}); 


		sdsl::sd_vector_builder builder(bwt_length,last_first.size());
		for(auto& idx: last_first){ builder.set(idx.first); }
		this->last = b_v(builder);

		std::cout << "Last built!" << std::endl;

		this->stpd = stpdA_; // set stpd array pointer
		first = sdsl::int_vector<>(last_first.size(),0,bitsize(uint64_t(bwt_length)));
		last_to_first = sdsl::int_vector<>(last_first.size(),0,bitsize(uint64_t(stpd->size())));
		std::cout << "Number of BWT runs indexed = " << last_first.size() << std::endl;


		{
			sdsl::int_vector<1> mapping(bwt_length,0);
			std::vector<uint_t> mapping_(stpd->size(),0);

			usafe_t i=0;
			for(const auto& e:*stpd){ mapping[e] = true; }
			sdsl::rank_support_v<1> rank1_(&mapping);
			for(const auto& e:*stpd){ mapping_[rank1_(e)] = i++; }

			i=0;
			usafe_t j=0, j_=0;
			f_bv = sdsl::int_vector<1>(last_first.size(),0);
			for(auto& entry: last_first)
			{
				if(mapping[entry.second+1])
				{
					f_bv[i] = true;
					last_to_first[j++] = mapping_[rank1_(entry.second+1)];
				}
				else { first[j_++] = entry.second; }
				i++;
			}
			this->first.resize(j_);
			this->last_to_first.resize(j);
			std::cout << "Number of SA samples in the STPD array = " << j << std::endl;
		}
	}

	uint_t phi_safe(const uint_t idx) const
	{
		if(idx != L)
		{
			auto res = last.successor_rank(idx);

			if(f_bv[res.second-1])
			{
				usafe_t no_ones = rank1_(res.second-1) + 1;
				return ((*stpd)[last_to_first[no_ones-1]]-1) - (res.first - idx);
			}
			else
			{
				usafe_t no_zeroes = res.second - rank1_(res.second);
				return first[no_zeroes-1] - (res.first - idx);
			}
		}
		else{ return -1; }
	}

	uint_t phi_unsafe(const uint_t idx) const
	{
		auto res = last.successor_rank(idx);

		if(f_bv[res.second-1])
		{
			usafe_t no_ones = rank1_(res.second-1) + 1;
			return ((*stpd)[last_to_first[no_ones-1]]-1) - (res.first - idx);
		}
		else
		{
			usafe_t no_zeroes = res.second - rank1_(res.second);
			return first[no_zeroes-1] - (res.first - idx);
		}
	}
	/*
	void test_phi(uint_t idx)
	{
		last.construct_rank_ds();
		last.construct_select_ds();
		rank1_ = sdsl::rank_support_v<1>(&f_bv);

		std::cout << idx << std::endl;
		while(idx != this->L)
		{
			idx = phi_safe(idx);
			std::cout << idx << std::endl;
		}
	}
	*/
	void load(std::istream& in)
	{
		in.read((char*)&L, sizeof(L));
		last.load(in);
		first.load(in);
	}

	uint_t serialize(std::ostream& out) const
	{
		uint_t w_bytes = 0;

		out.write((char*)&L, sizeof(L));

		w_bytes += sizeof(L);

		w_bytes += last.serialize(out);
		w_bytes += first.serialize(out);
		w_bytes += last_to_first.serialize(out);
		w_bytes += f_bv.serialize(out);

		w_bytes += sizeof(stpd);

		return w_bytes;
	}

private:
	// last samples
	b_v last;
	// last SA entry
	uint_t L;
	// first samples
	sdsl::int_vector<> first;
	sdsl::int_vector<> last_to_first;
	// first bitvector
	sdsl::int_vector<1> f_bv;
	sdsl::rank_support_v<1> rank1_;
	// pointer to stpa-array
	const sdsl::int_vector<>* stpd;

};

}

#endif // STPD_INDEX_PHI_INV_HPP_