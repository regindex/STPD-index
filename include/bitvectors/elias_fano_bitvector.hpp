/*
 * Construction of the Elias-Fano compressed bitvectors
 * 
 * This code is adapted from https://github.com/davidecenzato/extended_r-index/blob/main/sd_vector.hpp
 *
 */

#ifndef ELIAS_FANO_DS_HPP_
#define ELIAS_FANO_DS_HPP_

#include <cassert>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp>
#include <common.hpp>

namespace bv{

class elias_fano_bitvector{

public:
	// empty constructor
	elias_fano_bitvector(){}
	// parametrized constructor when input is a bitvector
	elias_fano_bitvector(std::vector<bool> &b)
	{
		if(b.size()==0) return;

		n = b.size();
		m = 0;

		sdsl::bit_vector bv_(b.size());

		for(uint64_t i=0;i<b.size();++i)
		{
			bv_[i] = b[i];
			if(b[i]){ m++; }
		}

		bv = sdsl::sd_vector<>(bv_);
	}
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
	// parametrized constructor when input is the vector builder
	elias_fano_bitvector(sdsl::sd_vector_builder& builder)
	{
		m = builder.items();
		// construct the compressed bitvector
		bv = sdsl::sd_vector<>(builder);
		// set bitvector len
		n = bv.size();
	}
	// parametrized constructor when input is filename containing the onset vector
	elias_fano_bitvector(const std::string& onset_file, const uint_t bsize, const int isize = 4)
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
        if (this != &other)
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

	// return the successor of i and its rank
	std::pair<uint_t,uint_t> successor_rank(uint_t i) const
	{
		std::pair<uint_t,uint_t> res;
		res.second = rank1_(i)+1;
		res.first = select1_(res.second);

		return res;
	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in)
	{
		in.read((char*)&n, sizeof(n));
		in.read((char*)&m, sizeof(m));
		bv.load(in);

		construct_rank_ds();
		construct_select_ds();
	}

	/* serialize the structure to the ostream
	 * \param out	 the ostream
	 */
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

private:
	// bitvector length
	uint_t n = 0;
	// no. bits set
	uint_t m = 0;
	// compressed bit vector
  	sdsl::sd_vector<> bv;
  	// rank data structure
  	sdsl::rank_support_sd<> rank1_;
  	// select data structure
  	sdsl::select_support_sd<> select1_;
};

}

#endif // ELIAS_FANO_DS_HPP_