#ifndef R_INDEX_PHI_HPP_
#define R_INDEX_PHI_HPP_

#include <elias_fano_bitvector.hpp>
#include <sdsl/int_vector.hpp>

namespace stpd{

template<class b_v = bv::elias_fano_bitvector>
class r_index_phi_inv
{
public:
	// empty constructor
	r_index_phi_inv(){}

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
				//std::cout << "c= " << curr << " " << curr_sa-1 << std::endl;
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

		sdsl::sd_vector_builder builder(bwt_length,last_first.size());
		for(auto& idx: last_first){ builder.set(idx.first); }
		last = b_v(builder);

		first = sdsl::int_vector<>(last_first.size(),0,bitsize(uint64_t(bwt_length)));
		if(verbose)
			std::cout << "Number of BWT runs indexed = " << first.size()+1 << std::endl;

		i=0;
		for(auto& entry: last_first)
		{
			first[i++] = entry.second;
			//std::cout << "(" << entry.first << "," << entry.second << ") " << std::endl;
		}
		//std::cout << "L: " << L << std::endl;

		bwt.close();
		sa.close();
	}

	int_t phi_safe(const uint_t idx) const
	{
		if(idx != L)
		{
			auto res = last.successor_rank(idx);
			//std::cout << res.first << " - " << res.second << std::endl;
			return first[res.second-1] - (res.first - idx);
		}
		else{ return -1; }
	}

	int_t phi_unsafe(const uint_t idx) const
	{
		auto res = last.successor_rank(idx);

		return first[res.second-1] - (res.first - idx);
	}

	void test_phi(uint_t idx)
	{
		last.construct_rank_ds();
		last.construct_select_ds();

		std::cout << idx << std::endl;
		while(idx != this->L)
		{
			idx = phi_safe(idx);
			std::cout << idx << std::endl;
		}
	}


	/*
	// parametrized constructor when input is the onset vector
	elias_fano_bitvector(const std::vector<uint_t>& onset, const uint_t bsize)
	{
		m = onset.size();
		// construct the compressed bitvector
		sdsl::sd_vector_builder builder(bsize,onset.size());
		for(auto idx: onset){ builder.set(idx); }
		bv = sdsl::sd_vector<>(builder);
		// set bitvector len
		n = bv.size();
	}
	// parametrized constructor when input is filename containing the onset vector
	elias_fano_bitvector(const std::string& onset, const uint_t bsize, const int isize = 4)
	{
		// open the input stream
		std::ifstream onset(onset_file,std::ifstream::binary);

		// compute onset size
		onset.seekg(0, std::ios::end);
		size_t onset_size = m = onset.tellg() / isize;
    	onset.seekg(0, std::ios::beg);
		// construct the compressed bitvector
		sdsl::sd_vector_builder builder(bsize,onset_size);
		// compute bitvector builder
		uint64_t currBitPos = 0;
		for(uint_t i=0;i<onset_size;++i){
			// get new onset pos
			onset.read(reinterpret_cast<char*>(&currBitPos), isize);
			builder.set(currBitPos);
		}
		bv = sdsl::sd_vector<>(builder);
		// close stream
		onset.close(); 
		// set bitvector len
		n = bv.size();

		// close stream 
		onset.close();
	}
	// copy constructor
	elias_fano_bitvector(const elias_fano_bitvector& other)
	{
		*this = other;
	}
	// overloading = operator
	elias_fano_bitvector& operator=(const elias_fano_bitvector & other)
	{
		if(this != &other)
		{
			n = other.n;
			m = other.m;
			bv = other.bv;
			rank1_ = other.rank1_;
			select1_ = other.select1_;
		}

	    return *this;
	}
	// move constructor
	elias_fano_bitvector(elias_fano_bitvector&& other) noexcept
	{
		*this = std::move(other);
	}
	// overloading = move operator
    elias_fano_bitvector& operator=(elias_fano_bitvector&& other)
    {
        if (this != &v)
        {
			n = other.n;
			m = other.m;
			bv = std::move(other.bv);
			rank1_ = std::move(other.rank1_);
			select1_ = std::move(other.select1_);
        }

        return *this;
    }

    // initialize rank data structure
	void construct_rank_ds()
	{
		assert(bv.size() > 0);
		sdsl::util::init_support(rank1_,&bv);
	}

	// initialize select data structure
	void construct_select_ds()
	{
		assert(bv.size() > 0);
		sdsl::util::init_support(select1_,&bv);
	}

	uint_t size() const { return n; }

	uint_t no_ones() const { return m; }
	
	uint_t rank1(uint_t i) const { return rank1_(i); }

	uint_t select1(uint_t i) const { return select1_(i+1); }

	uint_t at(uint_t i) const { return bv[i]; }

	uint_t gapAt(uint_t i) const
	{
		if(i==0){ return select1(0)+1; }
		return select1(i)-select1(i-1);
	}

	void load(std::istream& in)
	{
		in.read((char*)&n, sizeof(n));
		in.read((char*)&m, sizeof(m));
		bv.load(in);

		construct_rank_ds();
		construct_select_ds();
	}

	uint_t serialize(std::ostream& out) const
	{
		uint_t w_bytes = 0;

		out.write((char*)&n, sizeof(n));
		out.write((char*)&m, sizeof(m));

		w_bytes += sizeof(n) + sizeof(m);

		if(n == 0) return w_bytes;

		w_bytes += bv.serialize(out);

		return w_bytes;
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

		return w_bytes;
	}

private:
	// last samples
	b_v last;
	// last SA entry
	uint_t L;
	// first samples
	sdsl::int_vector<> first;

};

}

#endif // R_INDEX_PHI_HPP_