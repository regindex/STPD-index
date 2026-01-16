// Copyright (c) 2025, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 *  tabulated_binary_search_DNA: DNA optimized binary search with tabulation
 */

#pragma once

#include <cmath>
#include <common.hpp>

#include <elias_fano_static_sorted_map.hpp>

namespace stpd {

template< class text_oracle_ds = RLZ_DNA_sux<>,
          class elias_fano_ds  = stpd::ef::EliasFanoStaticSortedMap<> >
class tabulated_binary_search_DNA {

private:

	static constexpr uint8_t alph_w = 2; // DNA alphabet size width

	text_oracle_ds* oracle;  // random access text oracle (don't own this object)
	elias_fano_ds   samples; // Elias-Fano data structure

	int_t log_n, log_l; // SA and LCS sample widths

	uint_t  S;   // number of samples
	usafe_t N;   // text size
	int_t   len; // length of substrings stored in the Elias-Fano ds

public:

	void build( const std::string textFile, const std::string stpdArray, 
	            const std::string lcsArray, const std::string paArray,
                text_oracle_ds* oracle_, safe_t len_ = 15 )
	{
		// check tabulation length
		if(len_ < 1 or len_ > 31){
			std::cerr << "Tabulation size must be in the range [1..31]. exiting..."
					  << std::endl;
			exit(1);
		}

		// set private variables
		this->oracle = oracle_;
		this->N = oracle_->text_length();
		this->len = len_;

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

	    	std::cout << "		- STPD samples width = " << log_n << " bits per entry" << std::endl
	                  << "		- LCS values width = "   << log_l << " bits per entry" << std::endl;
		}

		// compute all STPD sample - lcs values pairs
		std::vector<std::pair<usafe_t,usafe_t>> key_value_pairs; 
		key_value_pairs.resize(S);
	    usafe_t a = 0, b = 0, c = 0, i = 0; 
	    safe_t min_lcs = std::numeric_limits<safe_t>::max();
	    std::vector<uchar_t> buffer(this->len,0);
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
					uchar_t c = oracle->extract(a);

					if(b > this->len){ b = this->len; }

					uint64_t samLcs = ((0ULL | a) << log_l) | b;
					key_value_pairs[i++].second = samLcs;
				}
				min_lcs = std::numeric_limits<safe_t>::max();
			}
		}
		file_stpd.close();
		file_lcs.close();
		file_pa.close();

		std::cout << "		- STPD array size = " << this->S <<
		std::endl << "		- Text size = " << this->N <<
		std::endl << "		- N/S = " << double(this->N)/S << std::endl;

		{ // Construct Elias-Fano binary search data structure
			usafe_t offset = 0, i = 0;
			std::ifstream file_text(textFile, std::ios::binary);
			if (!file_text.is_open()){ std::cerr << "Error: Could not open " << textFile << std::endl; exit(1); }

			safe_t curr = 0;
			for(i=0; i<S; ++i)
			{
				curr = key_value_pairs[i].second >> log_l;

				std::string text_buffer(this->len,'A');
				safe_t beg = std::max(static_cast<safe_t>(0),curr-this->len+1);
				safe_t len_s = std::min(static_cast<safe_t>(this->len),curr+1);

				file_text.seekg(beg, std::ios::beg);
				file_text.read(&text_buffer[this->len-len_s], len_s);
			  	file_text.clear();

			  	offset = bitpack_DNA_string(text_buffer,0,text_buffer.size());
				key_value_pairs[i].first = offset;
				offset = 0;
			}

			samples.build(key_value_pairs);

			file_text.close();
		}
	}
	
	// match all prefixes up to this->len
	std::pair<usafe_t,safe_t>
	match_first_prefix_with_tabulation(const std::string& pattern) const
	{
		usafe_t m = pattern.size(), 
        		i = std::min(static_cast<usafe_t>(this->len), m),
        		to_match = i;
		safe_t occ = -1;
		usafe_t bp_p = bitpack_DNA_string(pattern,0,i);

		while(i > 0) 
		{
			auto j = match_prefix(bp_p,i);
			occ = std::get<0>(j);

			bool is_leftmost = 
			std::min(static_cast<usafe_t>(this->len), i) > std::get<1>(j);
			if((occ != -1) and is_leftmost)
			{
				if(i < m)
				{
					usafe_t f = oracle->LCP(pattern,i,occ+1);

					i = i + f + 1;
					occ = occ + f;
				}
				else{ i++; }
				break;
			}

			i -= 1;
			bp_p = (bp_p & (1ULL << (this->len - 1) * alph_w) - 1) << alph_w;
		}

		// extend until we match at least this->len characters or
		// consume all the pattern
		while(i-1 < to_match)
		{
			auto j = match_prefix(pattern,0,i);
			occ = std::get<0>(j);

			usafe_t f = oracle->LCP(pattern,i,occ+1);
			i = i + f + 1;
			occ = occ + f;
		}

		return std::make_pair(i,occ);
	}

	// match all prefixes longer than this->len (approximate)
	safe_t 
	binary_search_with_tabluation(const std::string& P, usafe_t b, usafe_t e) const
	{
		usafe_t plen = e - b;
		assert(plen >= this->len);
		usafe_t to_match = this->len;

		// bitpack first to_match characters
		usafe_t search = bitpack_DNA_string(P, e-to_match, to_match);
		// search a range in the stpd array based on the fitst len characters suffix
		auto res = find_first_matching_sample_with_range(search);
		
		if(std::get<2>(res) < 0){ return -1; }

		safe_t sample = std::get<2>(res) >> log_l, lcs = to_match;
		usafe_t low = std::get<0>(res), high = std::get<1>(res);

		if( high > low+1 )
		{
			high -= 1;
			usafe_t mid = (low + high)/2;

			while( high > low )
			{	
				auto j = oracle->LCS_char(P, e-1, samples.get_value(mid) >> log_l); 
				
				if(j.first == (e - b)){ high = mid; }
				else if(j.second > P[e-j.first-1]) { high = mid-1; }
				else { low = mid+1; }

				mid = (low+high)/2;
			}

			sample = samples.get_value(low) >> log_l;
		}

		return sample;
	}

	bool is_bad_anchor(const std::string& pattern,usafe_t i,usafe_t anchor) const
	{
		return (oracle->LCS(pattern,i-1,anchor) != i);
	}

	usafe_t serialize(std::ostream& out)
	{
		usafe_t w_bytes = 0;

		out.write((char*)&N, sizeof(N));
		out.write((char*)&S, sizeof(S));
		out.write((char*)&len, sizeof(len));
		out.write((char*)&log_n, sizeof(log_n));
		out.write((char*)&log_l, sizeof(log_l));
		w_bytes += sizeof(N) + sizeof(S) + sizeof(len) +
		           		       sizeof(log_n) + sizeof(log_l);

		w_bytes += samples.serialize(out);

		return w_bytes;
	}

	void load(std::istream& in, text_oracle_ds* oracle_)
	{
		in.read((char*)&N, sizeof(N));
		in.read((char*)&S, sizeof(S));
		in.read((char*)&len, sizeof(len));
		in.read((char*)&log_n, sizeof(log_n));
		in.read((char*)&log_l, sizeof(log_l));

		oracle = oracle_;
		samples.load(in);
	}

private:

	int64_t find_first_matching_sample(uint64_t key, uint8_t key_width) const
	{
		auto r = samples.successor(key); 
		uint64_t s = r.first ^ key;

		uint8_t mbits = (s == 0) ? key_width : 
								   (__builtin_clzll(s) & ~1) - (64 - samples.universe_width());

		return (key_width <= mbits) ? r.second : -1;
	}

	std::tuple<usafe_t,usafe_t,safe_t>
	find_first_matching_sample_with_range(uint64_t key) const
	{
		auto s = samples.key_range(key); 
		
		if(std::get<0>(s) == std::get<1>(s)) return std::make_tuple(0,0,-1);
		
		if(std::get<2>(s) < this->len - 1)
		{
			uint64_t r = std::get<0>(s) + 1;
			auto p = samples.select(r);

			if(p.first != key) return std::make_tuple(0,0,-1);

			uint64_t r_ = std::get<1>(s);
			if(r == r_){ r_++; }

			return std::make_tuple(r,r_,p.second);
		}

		return s;
	}

	int64_t find_ith_matching_sample(uint64_t key, uint8_t key_width, uint64_t i) const
	{
		uint64_t r = samples.rank(key) + (i-1);
		auto p = samples.select(r);
		uint64_t s = p.first ^ key;

		uint8_t mbits = (s == 0) ? key_width : 
								   (__builtin_clzll(s) & ~1) - 
								   (64 - samples.universe_width());

		return (key_width <= mbits) ? p.second : -1;
	}

	usafe_t bitpack_DNA_string(const std::string& p, usafe_t beg, uchar_t plen) const
	{
		assert(this->len >= plen);

	    usafe_t res = 0;
	    for(uint_t i=0; i<plen; ++i)
	    {
	    	uchar_t c = p[beg+plen-i-1];
	        res = (res << alph_w) | 
	              static_cast<usafe_t>(dna_to_code_table[c]);
	    }
	    res = (res << (this->len - plen) * alph_w);

	    return res;
	}

	std::tuple<safe_t,safe_t> // return sample and lcs
	match_prefix(const std::string& P, usafe_t b, usafe_t e) const
	{
		usafe_t p_len = e - b;
		usafe_t to_match = std::min(static_cast<usafe_t>(this->len), p_len);

		// bitpack first to_match characters
		usafe_t search = bitpack_DNA_string(P, e-to_match, to_match);
		// search for the bitpacked pattern
		safe_t val = find_first_matching_sample(search, to_match * alph_w);

		// return an empty lower bound if we didn't match the pattern
		if(val < 0) return std::make_tuple(-1,-1);	

		usafe_t sample = val >> log_l, lcs = val & ((1ULL << log_l) - 1);
		// handle case where the sample is smaller than this->len
		if(sample+1 < to_match)
		{
			val = find_ith_matching_sample(search, to_match * alph_w, 2);

			if(val < 0) return std::make_tuple(-1,-1);

			sample = val >> log_l;
			lcs = val & ((1ULL << log_l) - 1);
		}

		return std::make_tuple(sample,lcs);
	}

	std::tuple<safe_t,safe_t> // return sample and lcs
	match_prefix(uint64_t p, usafe_t p_len) const
	{
		usafe_t to_match = std::min(static_cast<usafe_t>(this->len), p_len);

		// search for the bitpacked pattern
		safe_t val = find_first_matching_sample(p, to_match * alph_w);

		// return an empty lower bound if we didn't match the pattern
		if(val < 0) return std::make_tuple(-1,-1);	

		usafe_t sample = val >> log_l, lcs = val & ((1ULL << log_l) - 1);
		// handle case where the sample is smaller than this->len
		if(sample+1 < to_match)
		{
			val = find_ith_matching_sample(p, to_match * alph_w, 2);

			if(val < 0) return std::make_tuple(-1,-1);

			sample = val >> log_l;
			lcs = val & ((1ULL << log_l) - 1);
		}

		return std::make_tuple(sample,lcs);
	} 
};

}