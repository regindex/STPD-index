// Copyright (c) 2025, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 *  stpd-index: implementation of the Suffix Tree Path Decomposition index
 */

#ifndef STPD_INDEX_HPP_
#define STPD_INDEX_HPP_

#include <chrono>
#include <future>
#include <malloc_count.h>
// inverse phi functions
#include <r-index_phi_inv.hpp>
#include <stpd-index_phi_inv.hpp>
// bitpacked text oracle
#include <bitpacked_text_oracle.hpp>
// rlz text orcale
#include <RLZ_DNA.hpp>
// stpd-array binary search
#include <stpd_array_binary_search.hpp>
// stpd-array binary search opt
#include <stpd_array_binary_search_opt.hpp>

namespace stpd{

template<class STPDArray, class textOracle, class phiFunction>
class stpd_index{

private:
	// phi function data structure
	phiFunction phi;
	// text oracle data structure
	textOracle O;
	// stpd array search data structure
	STPDArray S;

public:
	// empty constructor
	stpd_index(){}

	void build_colex_m(const std::string &text_filepath, const std::string &sampling_filepath,
		               const std::string &rbwt_filepath, const std::string &pa_filepath)
	{
		std::cout << "Constructing the STPD-index for " << text_filepath << std::endl;
		std::cout << "Step 1) Constructing the random-access text oracle..." << std::endl;
		O.build(text_filepath);
		std::cout << "Step 2) Constructing the STPD-array binary search data structure..." << std::endl;
		S.build(sampling_filepath,&O,false); 
		std::cout << "Step 3) Constructing the phi function..." << std::endl;
	  	phi.build(rbwt_filepath,pa_filepath);
	  	
	  	std::cout << "Done!" << std::endl;
	}

	void build_colex_m(const std::string &text_filepath, const std::string &sampling_filepath,
		               const std::string &rbwt_filepath, const std::string &pa_filepath,
		               const std::string &lcs_filepath)
	{
		std::cout << "Constructing the STPD-index for " << text_filepath << std::endl;
		std::cout << "Step 1) Constructing the random-access text oracle..." << std::endl;
		O.build(text_filepath);
		std::cout << "Step 2) Constructing the STPD-array binary search data structure..." << std::endl;
		S.build(text_filepath,sampling_filepath,lcs_filepath,pa_filepath,&O,false); 
		std::cout << "Step 3) Constructing the phi function..." << std::endl;
	  	phi.build(rbwt_filepath,pa_filepath);
	  	
	  	std::cout << "Done!" << std::endl;
	}

	/*
	void build_colex_m_v2(const std::string &text_filepath)
	{
		std::cout << "Constructing the STPD-index for " << text_filepath << std::endl;
		std::cout << "Step 1) Constructing the random-access text oracle..." << std::endl;
		O.build(text_filepath);
		std::cout << "Step 2) Constructing the STPD-array binary search data structure..." << std::endl;
		S.build(text_filepath+".colex_m",&O); 
		std::cout << "Step 3) Constructing the phi function..." << std::endl;
	  	phi.build(text_filepath+".rbwt",text_filepath+".pa",S.get_array_point(),S.get_PA_size());
	} */

	// build suffixient index by indexing the supermaximal extensions
	void build_colex_pm(const std::string &text_filepath, const std::string &sampling_filepath,
		                const std::string &rbwt_filepath, const std::string &pa_filepath)
	{
		std::cout << "Constructing the STPD-index for " << text_filepath << std::endl;
		std::cout << "Step 1) Constructing the random-access text oracle..." << std::endl;
		O.build(text_filepath);
		std::cout << "Step 2) Constructing the STPD-array binary search data structure..." << std::endl;
		S.build(sampling_filepath,&O,true); 
		std::cout << "Step 3) Constructing the phi function..." << std::endl;
	  	phi.build(rbwt_filepath,pa_filepath);

	  	std::cout << "Done!" << std::endl;
	}
	
	usafe_t store(const std::string &index_filepath)
	{
		std::ofstream out(index_filepath);

		usafe_t O_bytes = O.serialize(out);
		std::cout << "Storing random-access data structure, size = " << O_bytes << " bytes" << std::endl;
		usafe_t S_bytes = S.serialize(out);
		std::cout << "Storing STPD-array data structure, size = " << S_bytes << " bytes" << std::endl;
		usafe_t phi_size = phi.serialize(out);
		std::cout << "Storing phi data structure, size = " << phi_size << " bytes" << std::endl;

		std::cout << "Index succesfully stored!" << std::endl;
		
		out.close();

		return O_bytes + S_bytes + phi_size;
	}
	
	void load(const std::string &index_filepath)
	{
		std::ifstream in(index_filepath);

		std::cout << "Loading the random-access text oracle..." << std::endl;
		O.load(in);
		std::cout << "Loading the STPD-array data structure..." << std::endl;
		S.load(in,&(this->O));
		std::cout << "Loading the phi data structure..." << std::endl;
		phi.load(in);

		std::cout << "Loading done!" << std::endl;

		in.close();
	}

	std::pair<std::vector<usafe_t>,double> 
						 locate_pattern_exp_search(const std::string &pattern) const
	{
		auto start = std::chrono::high_resolution_clock::now();

		usafe_t m = pattern.size();
		auto i_occ = this->S.locate_first_prefix(pattern);
		bool_t mismatch_found;

		while(i_occ.first-1 < m)
		{
			auto j = this->S.binary_search_lower_bound(pattern,0,i_occ.first);
			mismatch_found = std::get<2>(j);

			if(mismatch_found)
				return std::make_pair(std::vector<usafe_t>{},0);

			i_occ.second = std::get<0>(j);
			usafe_t f = O.LCP(pattern,i_occ.first,i_occ.second+1);
			i_occ.first = i_occ.first + f + 1;
			i_occ.second = i_occ.second + f;
		}

		usafe_t high = 2, low = 0;
		std::vector<usafe_t> res{usafe_t(i_occ.second)};
		while(true)
		{
			usafe_t phi_steps = high/2;
			while(phi_steps-- > 0)
			{
				i_occ.second = phi.phi_safe(i_occ.second);  // 0 1 2 3    
				if(i_occ.second == -1)
				{
					high -= phi_steps;
					binary_search(low,high,m,pattern,res);
					res.resize(low);

					std::chrono::duration<double> duration = 
							std::chrono::high_resolution_clock::now() - start;

					return std::make_pair(res,duration.count());			
				}

				res.push_back(i_occ.second);
			}

			usafe_t f = O.LCS(pattern,m-1,res[high-1]);
			if(f < m)
				break;

			low = high;
			high *= 2;
		}

		binary_search(low,high,m,pattern,res);
		res.resize(low);

		std::chrono::duration<double> duration = 
				std::chrono::high_resolution_clock::now() - start;

		return std::make_pair(res,duration.count());
	}

	std::pair<std::vector<usafe_t>,double> 
						 locate_pattern(const std::string pattern) const
	{
		auto start = std::chrono::high_resolution_clock::now();

		usafe_t i = 1, m = pattern.size();
		int_t lower_occ, upper_occ;
		bool_t mismatch_found;
		/*
		auto t1 = std::async(std::launch::async, [this, &pattern, &lower_occ, &mismatch_found]() {
		    this->lower_sample(pattern, lower_occ, mismatch_found);
		});
		auto t2 = std::async(std::launch::async, [this, &pattern, &upper_occ]() {
		    this->upper_sample(pattern, upper_occ);
		});
		t1.get(); t2.get();
		*/
		this->lower_sample(pattern, lower_occ, mismatch_found);
		this->upper_sample(pattern, upper_occ);

		if(mismatch_found)
			return std::make_pair(std::vector<usafe_t>{},0);

		std::vector<usafe_t> res{usafe_t(lower_occ)};
		while(lower_occ != upper_occ)
		{
			lower_occ = phi.phi_unsafe(lower_occ);
			res.push_back(lower_occ);
		}
		
		std::chrono::duration<double> duration = 
				std::chrono::high_resolution_clock::now() - start;

		return std::make_pair(res,duration.count());
	}

	//template<
    //std::pair<std::vector<usafe_t>, double> (*locate_all)(const std::string&) const
	//>
	void locate_fasta(const std::string patternFile) const
	{
		std::ifstream patterns(patternFile);
		std::ofstream   output(patternFile+".occs");

		std::string line, header;
		usafe_t i=0, c=0;
		std::pair<std::vector<usafe_t>,double> o;
		double tot_duration = 0;

		malloc_count_reset_peak();

		while(std::getline(patterns, line))
		{
			if(i%2 != 0)
			{

				if(this->S.is_index_large())
					o = locate_pattern(line);
				else
					o = locate_pattern_exp_search(line);

				output << header << std::endl;
				if(o.first.size() >= 0)
				{ 
					for(const auto& e:o.first)
						output << e << " ";
					output << std::endl;
				}
				else{ output << "-1 " << std::endl; }

				tot_duration += o.second;
				c += line.size();

			}
			else{ header = line; }
			i++;
		}

		patterns.close();
		output.close();

		std::cout << "Memory peak while running pattern matching queries = " <<
				     malloc_count_peak() << " bytes" << std::endl
		          << "Elapsed time while running pattern matching queries = " <<
				     tot_duration << " sec" << std::endl
		          << "Number of patterns = " << i/2 
		 		  << ", Total number of characters = " << c << std::endl
		          << "Elapsed time per pattern = " <<
				     (tot_duration/(i/2))*1000 << " milliSec" << std::endl
		          << "Elapsed time per character = " <<
				     (tot_duration/(c))*1000000 << " microSec" << std::endl;
	}

private:

	inline void binary_search(usafe_t& low, usafe_t& high, usafe_t m, 
		                      const std::string& pattern, const std::vector<usafe_t>& res) const
	{
		usafe_t mid = (low+high)/2;
		while( low < high )
		{		
			usafe_t f = O.LCS(pattern,m-1,res[mid]); 

			if(f == m)    
			{
				low = mid+1;
			}
			else
			{
				high = mid;
			}
 			
			mid = (low+high)/2;
		}
	}

	void lower_sample(const std::string& pattern, int_t& occ, bool_t& mismatch_found) const
	{
		usafe_t i = 1, m = pattern.size();
		occ = -1;

		while(i-1 < m)
		{
			auto j = this->S.binary_search_lower_bound(pattern,0,i);
			mismatch_found = std::get<2>(j);

			if(mismatch_found)
				return;

			occ = std::get<0>(j);
			usafe_t f = O.LCP(pattern,i,occ+1);
			i = i + f + 1;
			occ = occ + f;
		}
	}

	void upper_sample(const std::string& pattern, int_t& occ) const
	{
		usafe_t i = 1, m = pattern.size();
		occ = -1;

		while(i-1 < m)
		{
			auto j = this->S.binary_search_upper_bound(pattern,0,i);

			occ = j;
			usafe_t f = O.LCP(pattern,i,occ+1);
			i = i + f + 1;
			occ = occ + f;
		}
	}

	
};

}

#endif