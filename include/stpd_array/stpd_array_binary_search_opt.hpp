// Copyright (c) 2025, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 *  stpd_array_binary_search_opt: Optimized stpd array binary search implementation
 */

#ifndef STPD_ARRAY_BINARY_SEARCH_OPT_HPP_
#define STPD_ARRAY_BINARY_SEARCH_OPT_HPP_

#include <cmath>
#include <common.hpp>

#include <elias_fano_intlv.hpp> // elias fano dictionary data structure

namespace stpd{

template< class text_oracle_ds = RLZ_DNA_sux<>,
          class elias_fano_ds = sux::bits::InterleavedEliasFano<> >
class stpd_array_binary_search_opt {

public:

	stpd_array_binary_search_opt(){ }

	void build(const std::string textFile, const std::string stpdArray, 
		       const std::string lcsArray, const std::string paArray,
	           text_oracle_ds* O_,
	           bool_t large_ = false, safe_t len_ = 15, bool_t verbose = true)
	{
		{ // set input parameters
			this->large = large_;
			this->O = O_;
			this->N = O_->text_length();
			this->len = len_;
		}
		// compute the number of STPD samples
	   	std::ifstream file_stpd(stpdArray, std::ios::binary);
	   	if (!file_stpd.is_open()){ std::cerr << "Error: Could not open " << stpdArray << std::endl; exit(1); }
	   	file_stpd.seekg(0, std::ios::end);
	   	this->S = file_stpd.tellg()/5;
	   	file_stpd.seekg(0, std::ios::beg);
	   	// open PA and LCS files
	   	std::ifstream file_lcs(lcsArray, std::ios::binary);
	   	if (!file_lcs.is_open()){ std::cerr << "Error: Could not open " << lcsArray << std::endl; exit(1); }
	   	std::ifstream file_pa (paArray, std::ios::binary);
	   	if (!file_pa.is_open()){ std::cerr << "Error: Could not open " << paArray << std::endl; exit(1); }
	   	{ // compite width of samples and lcs entries
		    this->log_n = bitsize(this->N);
		    this->log_l = bitsize(this->len);

		    if(verbose)
		    	std::cout << "STPD samples width = " << log_n << " bits per entry" << std::endl
		                  << "LCS values width = " << log_l << " bits per entry" << std::endl;
		}

		// compute all STPD sample - lcs values pairs
		std::vector<std::pair<usafe_t,usafe_t>> key_value; 
		key_value.resize(S);
	    usafe_t a = 0, b = 0, c = 0, i = 0; 
	    safe_t min_lcs = std::numeric_limits<safe_t>::max();
	    std::vector<uchar_t> buffer(15,0);
		{ // read and bitpack the sample and lcs values
			while (file_stpd.read(reinterpret_cast<char*>(&buffer[0]), 5))
			{
				a = get_5bytes_uint(&buffer[0]);
				while (file_lcs.read(reinterpret_cast<char*>(&buffer[5]), 5) and
					    file_pa.read(reinterpret_cast<char*>(&buffer[10]), 5))
				{
					b = get_5bytes_uint(&buffer[5]);
					c = get_5bytes_uint(&buffer[10]);
					min_lcs = std::min(min_lcs,static_cast<safe_t>(b));

					if(c-1 == a)
						break;
				}
				if(a < this->N)
				{
					uchar_t c = O->extract(a);

					if(b > this->len){ b = this->len; }

					uint64_t samLcs = ((0ULL | a) << log_l) | b;
					key_value[i++].second = samLcs;
				}
				min_lcs = std::numeric_limits<safe_t>::max();
			}
		}
		file_stpd.close();
		file_lcs.close();
		file_pa.close();

		if(verbose) 
		{
			std::cout << "STPD array size = " << this->S <<
			std::endl << "Text size = " << this->N <<
			std::endl << "N/S = " << double(this->N)/S << std::endl;
		}

		{ // Construct Elias-Fano binary search data structure
			usafe_t offset = 0, i = 0;
			std::ifstream file_text(textFile, std::ios::binary);
			if (!file_text.is_open()){ std::cerr << "Error: Could not open " << textFile << std::endl; exit(1); }

			safe_t curr = 0;
			for(i=0; i<S; ++i)
			{
				curr = key_value[i].second >> log_l;

				std::string text_buffer(this->len,'A');
				safe_t beg = std::max(static_cast<safe_t>(0),curr-this->len+1);
				safe_t len_s = std::min(static_cast<safe_t>(this->len),curr+1);

				file_text.seekg(beg, std::ios::beg);
				file_text.read(&text_buffer[this->len-len_s], len_s);
			  	file_text.clear();

			  	bitpack_uint_DNA(reinterpret_cast<uint8_t*>(&offset),text_buffer);
				key_value[i].first = offset;
				offset = 0;
			}
			// compute the Elias-Fano data structure
			ef.build(key_value,pow(SIGMA_DNA,this->len),log_n+log_l);

			file_text.close();
		}
	}

	usafe_t sA_size() const { return this->S; }
	safe_t get_len() const { return this->len; }
	bool_t is_index_large() const { return this->large; }

	usafe_t serialize(std::ostream& out)
	{
		usafe_t w_bytes = 0;

		out.write((char*)&N, sizeof(N));
		out.write((char*)&S, sizeof(S));
		out.write((char*)&len, sizeof(len));
		out.write((char*)&large, sizeof(large));
		out.write((char*)&log_n, sizeof(log_n));
		out.write((char*)&log_l, sizeof(log_l));
		w_bytes += sizeof(N) + sizeof(S) + sizeof(len) + sizeof(large) +
		           						   sizeof(log_n) + sizeof(log_l);

		w_bytes += ef.serialize(out);

		return w_bytes;
	}

	void load(std::istream& in, text_oracle_ds* O_)
	{
		in.read((char*)&N, sizeof(N));
		in.read((char*)&S, sizeof(S));
		in.read((char*)&len, sizeof(len));
		in.read((char*)&large, sizeof(large));
		in.read((char*)&log_n, sizeof(log_n));
		in.read((char*)&log_l, sizeof(log_l));

		O = O_;
		ef.load(in);
	}

	// match all prefixes up to this->len
	std::pair<usafe_t,safe_t> locate_first_prefix(const std::string& pattern) const
	{
		usafe_t m = pattern.size(), 
        		i = std::min(static_cast<usafe_t>(this->len),m),
        		to_match = i;
		safe_t occ = -1;

		while(i > 0)
		{
			auto j = this->Elias_Fano_search_lower_bound(pattern,0,i);
			occ = std::get<0>(j);

			if((occ != -1) and (std::get<1>(j) > std::get<2>(j)))
			{
				if(i < m)
				{
					usafe_t f = O->LCP(pattern,i,occ+1);

					i = i + f + 1;
					occ = occ + f;
				}
				else{ i++; }
				break;
			}
			i--;
		}
		// extend until we match at least this->len characters or
		// consume all the pattern
		while(i-1 < to_match)
		{
			auto j = this->Elias_Fano_search_lower_bound(pattern,0,i);
			occ = std::get<0>(j);

			usafe_t f = O->LCP(pattern,i,occ+1);
			i = i + f + 1;
			occ = occ + f;
		}

		return std::make_pair(i,occ);
	}

	usafe_t 
	binary_search_upper_bound(const std::string& P, usafe_t pstart, usafe_t pend) const
	{
		std::cerr << "Not yet implemented!" << std::endl;
		exit(1);

		return 0;
	}

	// match all prefixes longer than this->len
	std::tuple<uint_t,uint_t,bool_t> 
		binary_search_lower_bound(const std::string& P, usafe_t b, usafe_t e) const
	{
		usafe_t plen = e - b;
		assert(plen >= this->len);
		uint_t to_match = this->len;

		// bitpack first to_match characters
		uint_t search = 0;
		this->bitpack_uint_DNA(reinterpret_cast<uint8_t*>(&search), P, this->len, 
			                                                e-to_match, to_match);
		// search a range in the stpd array based on the fitst len characters suffix
		auto res = 
		ef.lower_upper_bound_exact(search);

		if(std::get<0>(res) < 0){ return std::make_tuple(-1,0,1); }

		// run the binary search if needed
		if(std::get<0>(res)+1 == std::get<1>(res))
		{
			usafe_t lcp = O->LCS(P, e-1-this->len, (std::get<2>(res) >> log_l)-this->len) + this->len;

			return std::make_tuple((std::get<2>(res) >> log_l),lcp,lcp != plen);
		}
		else
		{
			auto bs_res =
			ef.binary_search_text_oracle(P, b, e, std::get<0>(res), std::get<1>(res)-1, log_l, O);

			return std::make_tuple(std::get<0>(bs_res),std::get<1>(bs_res),std::get<1>(bs_res) != plen);
		}
	}
	
	const sdsl::int_vector<>& get_array() const { return this->stpd; }
	const sdsl::int_vector<>* get_array_point() const { return &this->stpd; }
	usafe_t get_PA_size() const { return this->N; }
	usafe_t get_heur_len() const { return this->len; }

private:

	void inline bitpack_uint_DNA(uchar_t* t, const std::string& p, uchar_t len, 
	                                                  usafe_t beg, uchar_t plen) const
	{
	    assert(plen <= len);
	    uchar_t size = len; 
	    uchar_t offset = 4-(size%4);
	    for(uint_t i=0;i<plen;++i)
	        t[((size-i-1)/4)] |= 
	        bit_mask_table[(dna_to_code_table[p[beg+plen-i-1]]*4)+((offset+i)%4)];
	}

	void inline bitpack_uint_DNA(uchar_t* t, const std::string& p) const
	{
	    uchar_t size = p.size(); 
	    uchar_t offset = 4-(size%4); 
	    for(uint_t i=0;i<size;++i)
	        t[((size-i-1)/4)] |= 
	        bit_mask_table[(dna_to_code_table[p[size-i-1]]*4)+((offset+i)%4)];
	}

	std::tuple<safe_t,usafe_t,usafe_t> 
	Elias_Fano_search_lower_bound(const std::string& P, usafe_t b, usafe_t e) const
	{
		usafe_t plen = e - b;
		usafe_t to_match = std::min(static_cast<usafe_t>(this->len), plen);

		// bitpack first to_match characters
		usafe_t search = 0;
		bitpack_uint_DNA(reinterpret_cast<uint8_t*>(&search), P, this->len, e-to_match, to_match);
		// search for the bitpacked pattern
		safe_t val = ef.lower_bound(search,to_match * alph_w);

		// return an empty lower bound if we didn't match the pattern
		if(val < 0) return std::make_tuple(-1,0,0);	

		usafe_t sample = val >> log_l, lcs = val & ((1ULL << log_l) - 1);
		// handle case where the sample is smaller than this->len
		if(sample+1 < to_match)
		{
			val = ef.lower_bound_offset(search,to_match * alph_w,1);

			if(val < 0) return std::make_tuple(-1,0,0);

			sample = val >> log_l;
			lcs = val & ((1ULL << log_l) - 1);
		}

		return std::make_tuple(sample,to_match,lcs);
	}

	static constexpr uint8_t alph_w = 2; // DNA alphabet size width

	text_oracle_ds* O; // random access text oracle
	elias_fano_ds ef;  // Elias-Fano binary search data structure

	int_t log_n, log_l; // samples and lcs entries widths

	uint_t S;  // number of samples
	usafe_t N; // text size
	int_t len; // length of substrings stored in the Elias-Fano ds

	bool_t large; // contains both the colex- and colex+ samples
};

}

#endif // stpd_array_binary_search_opt_HPP