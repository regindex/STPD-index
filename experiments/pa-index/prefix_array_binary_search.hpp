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

	//std::tuple<uint_t,uint_t,bool_t> 
	uint_t locate_lower_bound(const std::string& pattern) const
	{
		// initialize binary search parameters
		uint_t low, mid, high, plen;
		//int_t lcp_low, lcp_high, lcp_mid;
		plen = pattern.size();
		low  = this->alph[pattern[plen-1]];
		high = this->alph[pattern[plen-1]+1];
		/*
		usafe_t count = 0;
		for(usafe_t i=0;i<PA.size();++i)
			if(PA[i]==0) count++;
		std::cout << "-* " << count << std::endl;
		*/
		// stop if first pattern character doesn't occur in the text
		if((high - low) > 0)
			{ 
				if(plen == 1){ return low; }
				high--;
				//lcp_low = lcp_high = -1; 
				mid = (low+high)/2;
			}
		else{ return -1; }


		while( high != low )
		{
			//std::cout << "low = "<< low << " high = " << high << std::endl;
			//std::cout << "mid = " << mid << " -> " << PA[mid] << std::endl;		
			auto j = O->LCS_char(pattern,plen-1,this->PA[mid]); 

			//std::cout << j.first << "," << j.second << std::endl;
	
			//if(j.first == plen){ return std::make_tuple(this->PA[mid],plen,false); }

			if(j.first == plen or j.second > pattern[plen-j.first-1])
			{  
				high = mid;
			}
			else{
				low = mid + 1;
			}
			mid = (low+high)/2;
		}

		//if(low != high){ low++; }

		//std::cout << "--> low = "<< low << " "<< PA[low] << " " << PA[low+1] << std::endl;
		//exit(1);

		return low;
	}

	uint_t locate_upper_bound(const std::string& pattern) const
	{
		// initialize binary search parameters
		uint_t low, mid, high, plen;
		//int_t lcp_low, lcp_high, lcp_mid;
		plen = pattern.size();
		low  = this->alph[pattern[plen-1]];
		high = this->alph[pattern[plen-1]+1];

		// stop if first pattern character doesn't occur in the text
		if((high - low) > 0)
			{ 
				if(plen == 1){ return high; }
				//lcp_low = lcp_high = -1; 
				mid = (low+high)/2;
			}
		else{ return -1; }

		while( high-low > 1 )
		{		
			//std::cout << "low = "<< low << " high = " << high << std::endl;
			//std::cout << "mid = " << mid << " -> " << PA[mid] << std::endl;	
			auto j = O->LCS_char(pattern,plen-1,this->PA[mid]); 
	
			//if(j.first == plen){ return std::make_tuple(this->PA[mid],plen,false); }

			if(j.second > pattern[plen-j.first-1])
			{  
				high = mid;
			}
			else
			{
				low = mid;
			}
			mid = (low+high)/2;
		}

		//if(low == high){ high++; }

		//std::cout << std::endl << "--> low = "<< low << " "<< PA[low-1] << " " << PA[low] << std::endl;
		//std::cout << "--> high = "<< high << " "<< PA[high-1] << " " << PA[high] << " " << PA[high+1] << std::endl;
		//exit(1);

		return high;
	}

	std::vector<uint_t> get_SA_range(uint_t b, uint_t e) const
	{
		std::vector<uint_t> range(e-b,0);
		for(uint_t i=0;i<range.size();++i)
		{
			range[i] = PA[b+i];
		}

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