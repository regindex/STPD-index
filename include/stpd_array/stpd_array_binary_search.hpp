// Copyright (c) 2025, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 *  stpd_array_binary_search: Baseline stpd-array index implementation
 */

#ifndef STPD_ARRAY_BINARY_SEARCH_HPP_
#define STPD_ARRAY_BINARY_SEARCH_HPP_

#include <cmath>
#include <common.hpp>
#include <sdsl/int_vector.hpp>

namespace stpd{

template<class text_oracle>
class stpd_array_binary_search{

public:

	stpd_array_binary_search(){ }

	void build(const std::string stpdArray, text_oracle* O_, 
		             bool_t large_, bool_t verbose = true)
	{
		{
			this->large = large_;
			this->O = O_;
			this->N = O_->text_length();
		}
	   	std::ifstream file_stpd(stpdArray, std::ios::binary);
	   	file_stpd.seekg(0, std::ios::end);
	   	this->S = file_stpd.tellg()/5;
	   	file_stpd.seekg(0, std::ios::beg);
	   	{
		    int_t log_n = bitsize(this->N);
		    this->stpd = sdsl::int_vector<>(S,0,log_n);
		    if(verbose)
		    	std::cout << "STPD array width = " << log_n << " bits per entry"
		     << std::endl;
		}
	    usafe_t a = 0, i = 0; 
	    std::vector<uchar_t> buffer(5,0);
		{
			this->alph = sdsl::int_vector<>(SIGMA,0,40);
			while (file_stpd.read(reinterpret_cast<char*>(&buffer[0]), 5))
			{
				a = get_5bytes_uint(&buffer[0]);
				if(a < this->N)
				{
					uchar_t c = O->extract(a);
					this->stpd[i++] = a;
					this->alph[c]++;
				}
			}
		}
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
		/*
		std::cout << "-->AAA test" << std::endl;
		std::string pattern = "AAA";
		binary_search_lower_bound(pattern,0,3);
		binary_search_upper_bound(pattern,0,3);
		exit(1);
		*/
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
		w_bytes += sizeof(N) + sizeof(S) + sizeof(len) + sizeof(large);

		w_bytes += stpd.serialize(out);
		w_bytes += alph.serialize(out);

		return w_bytes;
	}

	void load(std::istream& in, text_oracle* O_)
	{
		std::cout << "start loading" << std::endl;
		O = O_;
		in.read((char*)&N, sizeof(N));
		in.read((char*)&S, sizeof(S));
		in.read((char*)&len, sizeof(len));
		in.read((char*)&large, sizeof(large));

		stpd.load(in);
		alph.load(in);
	}

	std::tuple<uint_t,uint_t,bool_t> 
		binary_search_lower_bound(const std::string& pattern,uint_t pstart,uint_t pend) const
	{
		// initialize binary search parameters
		uint_t low, mid, high, lcp, plen;
		plen = pend - pstart;
		//std::cout << "plen= " << plen << std::endl;
		low  = this->alph[pattern[pend-1]];
		high = this->alph[pattern[pend-1]+1];

		//std::cout << low << " - " << high << std::endl;

		// stop if first pattern character doesn't occur in the text
		if((high - low) > 0)
		{ 
			if(plen == 1)
			{
				//std::cout << "return " << low << " - " << high-1 << std::endl;
				//exit(1);
				return std::make_tuple(this->stpd[low],1,false);
			}
			high--;
			lcp = O->LCS(pattern,pend-1,stpd[high]); 
			mid = (low+high)/2;
		}
		else
			return std::make_tuple(-1,0,true);

		//std::cout << low << " - " << high << std::endl;

		while( low < high )
		{		
			auto j = O->LCS_char(pattern,pend-1,this->stpd[mid]); 
	
			//if(j.first == plen)
			//	return std::make_tuple(this->stpd[mid],plen,false);    -- B AAAA , C AAAA

			if((j.first != plen) and (j.second < pattern[pend-j.first-1]))    
			{
				low = mid+1;
			}
			else
			{
				high = mid;
				lcp = j.first;
			}
 			
 			//if(lcp_low < j.first){ lcp_low = j.first; }
			mid = (low+high)/2;
		}

		//if(lcp_low  == -1){ lcp_low  = O->LCS_char(pattern,pend-1,stpd[low]).first; }

		//std::cout << "ritorna: " << low << " - " << lcp << std::endl;

		return std::make_tuple(this->stpd[low],lcp,(lcp != plen));
	}

	uint_t 
		binary_search_upper_bound(const std::string& pattern,uint_t pstart,uint_t pend) const
	{
		// initialize binary search parameters
		uint_t low, mid, high, plen;
		plen = pend - pstart;
		//std::cout << "plen= " << plen << std::endl;
		low  = this->alph[pattern[pend-1]];
		high = this->alph[pattern[pend-1]+1];

		//std::cout << low << " - " << high << std::endl;

		// stop if first pattern character doesn't occur in the text
		if((high - low) > 0)
		{ 
			if(plen == 1)
			{
				//std::cout << "return " << low << " - " << high << std::endl;
				//exit(1);
				return this->stpd[high-1];
			}
			//high--;
			mid = (low+high)/2;
		}
		else
			return -1;


		while( low < high )
		{		
			//std::cout << low << " - " << high << " : " << mid << std::endl;
			auto j = O->LCS_char(pattern,pend-1,this->stpd[mid]); 
	
			//if(j.first == plen)
			//	return std::make_tuple(this->stpd[mid],plen,false);    -- B AAAA , C AAAA
			
			if((j.first != plen) and (j.second > pattern[pend-j.first-1])) 
			{
				high = mid;
			}
			else
			{
				low = mid+1;
			}
			 			
 			//if(lcp_low < j.first){ lcp_low = j.first; }
			mid = (low+high)/2;
			//std::cout << low << " - " << high << " : " << mid << std::endl;
		}

		//if(lcp_low  == -1){ lcp_low  = O->LCS_char(pattern,pend-1,stpd[low]).first; }

		//std::cout << "ritorna: " << high  << std::endl;

		return this->stpd[high-1];
	}

	const sdsl::int_vector<>& get_array() const { return stpd; }
	const sdsl::int_vector<>* get_array_point() const { return &stpd; }
	usafe_t get_PA_size() const { return N; }

private:

	text_oracle* O;

	sdsl::int_vector<> stpd;
	sdsl::int_vector<> alph;

	uint_t S;
	usafe_t N;
	int_t len;

	bool_t large;
};

}

#endif // stpd_array_binary_search_HPP_