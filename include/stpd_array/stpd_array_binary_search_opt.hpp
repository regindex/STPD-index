// Copyright (c) 2025, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 *  stpd_array_binary_search_opt: Optimized stpd-array index implementation
 */

#ifndef STPD_ARRAY_BINARY_SEARCH_OPT_HPP_
#define STPD_ARRAY_BINARY_SEARCH_OPT_HPP_

#include <cmath>
#include <common.hpp>
#include <elias_fano_bitvector.hpp>
#include <succinct_bitvector.hpp>
#include <sdsl/int_vector.hpp>

namespace stpd{

template<class text_oracle, class elias_fano_ds = bv::elias_fano_bitvector,
         class bitvector = bv::succinct_bitvector>
class stpd_array_binary_search_opt{

public:

	stpd_array_binary_search_opt(){ }

	void build(const std::string textFile, const std::string stpdArray, 
		       const std::string lcsArray, const std::string paArray,
	           text_oracle* O_, bool_t large_, bool_t verbose = true)
	{
		{
			this->large = large_;
			this->O = O_;
			this->N = O_->text_length();
			this->len = 15;
		}
	   	std::ifstream file_stpd(stpdArray, std::ios::binary);
	   	file_stpd.seekg(0, std::ios::end);
	   	this->S = file_stpd.tellg()/5;
	   	file_stpd.seekg(0, std::ios::beg);
	   	std::ifstream file_lcs(lcsArray, std::ios::binary);
	   	std::ifstream  file_pa(paArray, std::ios::binary);
	   	sdsl::int_vector<> LCS;
	   	{
		    this->log_n = bitsize(this->N);
		    this->log_l = bitsize(this->len);
		    this->stpd = sdsl::int_vector<>(S*(log_n+log_l)+1,0,1);
		    LCS = sdsl::int_vector<>(S,0,log_n);

		    if(verbose)
		    	std::cout << "STPD array width = " << log_n + log_l << " bits per entry"
		     << std::endl;
		}
	    usafe_t a = 0, b = 0, c = 0, i = 0; 
	    safe_t min_lcs = std::numeric_limits<safe_t>::max();
	    std::vector<uchar_t> buffer(15,0);
		{
			this->alph = sdsl::int_vector<>(SIGMA,0,40);
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
					LCS[i] = min_lcs;
					set_sample_lcs(i++,a,b);
					this->alph[c]++;
				}
				min_lcs = std::numeric_limits<safe_t>::max();
			}
		}
		file_stpd.close();
		file_lcs.close();
		file_pa.close();
	    {
		    usafe_t sum = 0;
		    for(i=0;i<SIGMA;++i){ usafe_t tmp = this->alph[i]; this->alph[i] = sum; sum += tmp; }
	  	} 
		if(verbose) 
		{
			std::cout << "STPD array size = " << this->S <<
			std::endl << "Text size = " << this->N <<
			std::endl << "N/S = " << double(this->N)/S << std::endl;
		}
		{ // Construct Elias-Fano optimization
			std::vector<uint_t> onset;
			uint_t offset = 0;
			std::ifstream file_text(textFile, std::ios::binary);
			if (!file_text.is_open())
			{
			std::cerr << "Error: Could not open " << 
						 			 textFile << std::endl; exit(1);
			}
			safe_t prev = -1;
			std::vector<bool_t> tmp_bv(this->stpd.size()+1,0);

			for(i=0;i<S;++i)
			{
				if(LCS[i] < this->len)
				{
					tmp_bv[i] = 1; 
					if(prev != -1)
					{
						std::string text_buffer(this->len,'A');
						safe_t beg = std::max(static_cast<safe_t>(0),prev-this->len+1);
						safe_t len_ = std::min(static_cast<safe_t>(this->len),prev+1);

						file_text.seekg(beg, std::ios::beg);
						file_text.read(&text_buffer[this->len-len_], len_);
					  	file_text.clear();

					  	bitpack_uint_DNA(reinterpret_cast<uint8_t*>(&offset),text_buffer);
						onset.push_back(offset);
						offset = 0;
					}
					prev = this->get_sample(i);
				}
			}
			tmp_bv[S] = 1;

			std::string text_buffer(this->len,'A');
			safe_t beg = std::max(static_cast<safe_t>(0),prev-this->len+1);
			safe_t len_ = std::min(static_cast<safe_t>(this->len),prev+1);
			file_text.seekg(beg, std::ios::beg);
			file_text.read(&text_buffer[this->len-len_], len_);

			bitpack_uint_DNA(reinterpret_cast<uint8_t*>(&offset),text_buffer);
			onset.push_back(offset);
			this->bv = bitvector(tmp_bv);
			this->ef = elias_fano_ds(onset,pow(SIGMA_DNA,this->len));

			file_text.close();
		}
	}

	usafe_t sA_size() const { return this->S; }
	int_t get_len() const { return this->len; }
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

		w_bytes += stpd.serialize(out);
		w_bytes += ef.serialize(out);
		w_bytes += bv.serialize(out);
		w_bytes += alph.serialize(out);

		return w_bytes;
	}

	void load(std::istream& in, text_oracle* O_)
	{
		O = O_;
		in.read((char*)&N, sizeof(N));
		in.read((char*)&S, sizeof(S));
		in.read((char*)&len, sizeof(len));
		in.read((char*)&large, sizeof(large));
		in.read((char*)&log_n, sizeof(log_n));
		in.read((char*)&log_l, sizeof(log_l));

		stpd.load(in);
		ef.load(in);
		bv.load(in);
		alph.load(in);
	}

	std::pair<usafe_t,int_t> locate_first_prefix(const std::string& pattern) const
	{
		usafe_t m = pattern.size(), 
        i = std::min(static_cast<usafe_t>(this->len),m);
		int_t occ = -1;

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

		return std::make_pair(i,occ);
	}

	uint_t 
		binary_search_upper_bound(const std::string& pattern,uint_t pstart,uint_t pend) const
	{
		// initialize binary search parameters
		uint_t low, mid, high, plen;
		plen = pend - pstart;
		low  = this->alph[pattern[pend-1]];
		high = this->alph[pattern[pend-1]+1];

		// stop if first pattern character doesn't occur in the text
		if((high - low) > 0)
		{ 
			if(plen == 1)
			{
				return this->get_sample(high-1);
			}
			mid = (low+high)/2;
		}
		else
			return -1;


		while( low < high )
		{		
			auto j = O->LCS_char(pattern,pend-1,this->get_sample(mid)); 
			
			if((j.first != plen) and (j.second > pattern[pend-j.first-1])) 
			{
				high = mid;
			}
			else
			{
				low = mid+1;
			}
			 			
			mid = (low+high)/2;
		}

		return this->get_sample(high-1);
	}

	std::tuple<uint_t,uint_t,bool_t> 
		binary_search_lower_bound(const std::string& pattern,uint_t pstart,uint_t pend) const
	{
		uint_t plen = pend-pstart;
		uint_t to_match = std::min(static_cast<uint_t>(this->len),plen);
		uint_t high, low, mid, lcp;
		uint_t r, r_;

		{ // Run the Elias-Fano search
			// bitpack first to_match characters
			uint_t search = 0;
			this->bitpack_uint_DNA(reinterpret_cast<uint8_t*>(&search), pattern, this->len, 
				                                                   pend-to_match, to_match);
			// search for the bitpacked pattern
			       r = ef.rank1(search);    
			uint_t s = ef.select1(r);
			uint_t mlen = (__builtin_clz(search ^ s)-((sizeof(uint_t)*8)-(this->len*2)))/2;

			if( mlen < to_match )
				return std::make_tuple(-1,0,true);
		
			if(to_match == this->len)
			{ 
				r_ = r+1;
			}
			else
			{
				uint_t search_ = search | (1 << (this->len-to_match)*2)-1;
				this->bitpack_uint_DNA(reinterpret_cast<uint8_t*>(&search_), pattern, this->len, 
					                                                   pend-to_match, to_match);
				r_ = ef.rank1(search_)-1;
			}
			low  = bv.select(r); high = bv.select(r_);

			if(this->get_sample(low)+1 < this->len)
			{
				r++;
				s = ef.select1(r);
				mlen = (__builtin_clz(search ^ s)-((sizeof(uint_t)*8)-(this->len*2)))/2;

				if( mlen < to_match )
					return std::make_tuple(-1,0,0);

				low  = bv.select(r);
				if(r == r_){ r_++; high = bv.select(r_); }
			}
		}
		{ // Run the binary search
			// stop if first pattern character doesn't occur in the text
			if((high - low) > 0)
			{ 
				high--;
				lcp = O->LCS(pattern,pend-1,this->get_sample(high)); 
				mid = (low+high)/2;
			}
			else
				return std::make_tuple(-1,0,true);

			while( low < high )
			{	
				auto j = O->LCS_char(pattern,pend-1,this->get_sample(mid)); 
		
				if((j.first != plen) and (j.second < pattern[pend-j.first-1]))    
				{
					low = mid+1;
				}
				else
				{
					high = mid;
					lcp = j.first;
				}
	 			
				mid = (low+high)/2;
			}
		}

		return std::make_tuple(this->get_sample(low),lcp,(lcp != plen));
	}
	
	const sdsl::int_vector<>& get_array() const { return this->stpd; }
	const sdsl::int_vector<>* get_array_point() const { return &this->stpd; }
	usafe_t get_PA_size() const { return this->N; }
	usafe_t get_heur_len() const { return this->len; }

private:

	void inline bitpack_uint_DNA(uchar_t* t, const std::string& p, uchar_t len, 
	                                                  uint_t beg, uchar_t plen) const
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

	inline void set_sample_lcs(usafe_t i, usafe_t sample, usafe_t lcs)
	{ 
		this->stpd.set_int(i*(log_n+log_l), sample, log_n);
		this->stpd.set_int(i*(log_n+log_l)+log_n, lcs, log_l);
	}

	inline void get_sample_lcs(usafe_t i, usafe_t& sample, usafe_t& lcs) const
	{ // 0 0 0 0 0 0 0 0 0
		usafe_t e = this->stpd.get_int(i*(log_n+log_l), log_n+log_l);
		sample = e & ((1 << log_n)-1);
		lcs    = e >> log_n;
	}

	inline usafe_t get_sample(usafe_t i) const
	{
		return this->stpd.get_int(i*(log_n+log_l), log_n);
	}

	inline usafe_t get_lcs(usafe_t i) const
	{
		return this->stpd.get_int(i*(log_n+log_l)+log_n, log_l);
	}

	std::tuple<safe_t,uint_t,uint_t> 
		Elias_Fano_search_lower_bound(const std::string& pattern,uint_t pstart,uint_t pend) const
	{
		uint_t plen = pend-pstart;
		uint_t to_match = std::min(static_cast<uint_t>(this->len),plen);

		// bitpack first to_match characters
		uint_t search = 0;
		this->bitpack_uint_DNA(reinterpret_cast<uint8_t*>(&search), pattern, this->len, 
			                                                   pend-to_match, to_match);
		// search for the bitpacked pattern
		uint_t r = ef.rank1(search);
		uint_t s = ef.select1(r);
		uint_t mlen = (__builtin_clz(search ^ s)-((sizeof(uint_t)*8)-(this->len*2)))/2;

		if( mlen < to_match )
			return std::make_tuple(-1,0,0);	

		usafe_t sample, lcs;
		this->get_sample_lcs(bv.select(r),sample,lcs);
		if(sample+1 < this->len)
		{
			r++;
			s = ef.select1(r);
			mlen = (__builtin_clz(search ^ s)-((sizeof(uint_t)*8)-(this->len*2)))/2;

			if( mlen < to_match )
				return std::make_tuple(-1,0,0);

			this->get_sample_lcs(bv.select(r),sample,lcs);
		}

		//return std::make_tuple(bv.select(r),to_match,false);
		return std::make_tuple(sample,to_match,lcs);
	}

	text_oracle* O;

	sdsl::int_vector<> stpd;
	sdsl::int_vector<> alph;

	elias_fano_ds ef;
	bitvector bv;

	int_t log_n, log_l;

	uint_t S;
	usafe_t N;
	int_t len;

	bool_t large;
};

}

#endif // stpd_array_binary_search_HPP_