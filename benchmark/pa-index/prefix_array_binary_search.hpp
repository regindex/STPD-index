// Copyright (c) 2025, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 *  prefix_array_binary_search: index performing binary search on the full prefix array
 */

#ifndef PREFIX_ARRAY_BINARY_SEARCH_HPP_
#define PREFIX_ARRAY_BINARY_SEARCH_HPP_

#include <cmath>
#include <common.hpp>
#include <sdsl/int_vector.hpp>

namespace stpd{

template<class text_oracle = RLZ_DNA_sux<>>
class prefix_array_binary_search{

public:

	prefix_array_binary_search(){ this->alph = sdsl::int_vector<>(SIGMA,0,40); }

	int_t get_len(){ return this->len; }

	void build(std::string input_text, std::string input_pa, text_oracle* O_, bool_t verbose = true)
	{
		this->O = O_;
	  	{
		  	std::ifstream file_text(input_text, std::ios::binary);
		    file_text.seekg(0, std::ios::end);
		    this->N = this->S = file_text.tellg();
		    file_text.seekg(0, std::ios::beg);
		    uchar_t c;
			while (file_text.read(reinterpret_cast<char*>(&c), 1)){ this->alph[c]++; }
		    file_text.close();
		}
		usafe_t a = 0, i = 0;
	    {
		    usafe_t sum = 0;
		    for(i=0;i<SIGMA;++i){ usafe_t tmp = this->alph[i]; this->alph[i] = sum; sum += tmp; }
	  	}
	   	{
		    int_t log_n = bitsize(this->S);
		    this->PA = sdsl::int_vector<>(this->S,0,log_n);
		    if(verbose)
		    	std::cout << "Prefix array size = " << log_n << " bits per entry"
		     << std::endl;
		}
	   	std::ifstream file_pa(input_pa, std::ios::binary);
		if (!file_pa.is_open())
		{
	    std::cerr << "Error: Could not open " << 
	    			 			 input_pa << std::endl; exit(1);
		}
		{
			i = 0;
			std::vector<uchar_t> buffer(5,0);
			file_pa.seekg(5, std::ios::beg);
			while (file_pa.read(reinterpret_cast<char*>(&buffer[0]), 5))
			{
				a = get_5bytes_uint(&buffer[0])-1;
				this->PA[i++] = a;
			}
		}
		file_pa.close();
		if(verbose)
		{
			std::cout << "Prefix array size = " << this->S <<
			std::endl << "Text size = " << this->N <<
			std::endl << "N/S = " << double(this->N)/S << std::endl;
		}
	}

	usafe_t serialize(std::ostream& out)
	{
		usafe_t w_bytes = 0;

		out.write((char*)&N, sizeof(N));
		out.write((char*)&S, sizeof(S));
		w_bytes += sizeof(N) + sizeof(S);

		w_bytes += PA.serialize(out);
		w_bytes += alph.serialize(out);

		return w_bytes;
	}

	void load(std::istream& in, text_oracle* O_)
	{
		O = O_;
		in.read((char*)&N, sizeof(N));
		in.read((char*)&S, sizeof(S));

		PA.load(in);
		alph.load(in);
	}

	usafe_t range_lower_boundary(const std::string& pattern) const
	{
		// initialize binary search parameters
		usafe_t low, mid, high, plen;
		plen = pattern.size();
		low  = this->alph[ pattern[plen - 1] ];
		high = this->alph[ pattern[plen - 1] + 1 ];

		// stop if first pattern character doesn't occur in the text
		if( (high - low) > 0 )
		{ 
			if(plen == 1)
				return low;

			high--;
			mid = (low + high) / 2;
		}

		// run the binary search
		while(high != low)
		{	
			auto lcs = O->LCS_char(pattern, plen-1, this->PA[mid]); 

			if(lcs.first == plen or 
			   lcs.second > pattern[plen - lcs.first - 1])
				high = mid;
			else
				low = mid + 1;

			mid = (low + high) / 2;
		}

		return low;
	}

	std::pair<usafe_t,safe_t>
	range_lower_boundary_lcs(const std::string& pattern) const
	{
		// initialize binary search parameters
		usafe_t low, mid, high, plen;
		safe_t curr_lcs = -1;
		plen = pattern.size();
		low  = this->alph[ pattern[plen - 1] ];
		high = this->alph[ pattern[plen - 1] + 1 ];

		// stop if first pattern character doesn't occur in the text
		if( (high - low) > 0 )
		{ 
			if(plen == 1)
				return std::make_pair(low,0);

			high--;
			mid = (low + high) / 2;
		}

		// run the binary search
		while(high != low)
		{	
			auto lcs = O->LCS_char(pattern, plen-1, this->PA[mid]); 

			if(lcs.first == plen or 
			   lcs.second > pattern[plen - lcs.first - 1])
			{
				high = mid;
				curr_lcs = lcs.first;
			}
			else
				low = mid + 1;

			mid = (low + high) / 2;
		}

		if(curr_lcs < 0)
			curr_lcs = O->LCS_char(pattern, plen-1, this->PA[low]).first;

		return std::make_pair(low,curr_lcs);
	}

	usafe_t range_upper_boundary(const std::string& pattern) const
	{
		// initialize binary search parameters
		usafe_t low, mid, high, plen;
		plen = pattern.size();
		low  = this->alph[ pattern[plen - 1] ];
		high = this->alph[ pattern[plen - 1] + 1 ];

		// stop if first pattern character doesn't occur in the text
		if( (high - low) > 0 )
		{ 
			if(plen == 1)
				return high;

			mid = (low + high) / 2;
		}

		// run the binary search
		while(high-low > 1)
		{		
			auto lcs = O->LCS_char(pattern, plen-1, this->PA[mid]); 

			if(lcs.second > pattern[plen - lcs.first - 1])
				high = mid;
			else
				low = mid;

			mid = (low + high) / 2;
		}

		return high;
	}

	std::vector<uint_t> get_SA_range(uint_t b, uint_t e) const
	{
		std::vector<uint_t> range(e-b,0);
		for(usafe_t i=0;i<range.size();++i)
		{
			range[i] = PA[b+i];
		}

		return range;
	}

	std::vector<uint_t> get_SA_range_thr(uint_t b, uint_t e, usafe_t thr) const
	{
		std::vector<uint_t> range( std::min(static_cast<usafe_t>(e - b), thr), 0 );
		
		for(usafe_t i=0; i<range.size(); ++i)
			range[i] = PA[b+i];

		return range;
	}

	uint_t operator[](uint_t i) const
	{
		return PA[i];
	}

private:

	text_oracle* O;

	sdsl::int_vector<> PA;
	sdsl::int_vector<> alph;

	uint_t S;
	usafe_t N;
	int_t len;
};

}

#endif