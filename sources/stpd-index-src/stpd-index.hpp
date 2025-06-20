// Copyright (c) 2025, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 *  stpd-index: implementation of the Suffix Tree Path Decomposition index
 */

#ifndef STPD_INDEX_HPP_
#define STPD_INDEX_HPP_

#include <chrono>
#include <malloc_count.h> 

#include <r-index_phi_inv_intlv.hpp> // phi function
#include <RLZ_DNA_sux.hpp> // rlz random access text orcale
#include <stpd_array_binary_search.hpp> // binary search ds
#include <stpd_array_binary_search_opt.hpp> // optimized binary search ds

namespace stpd{

template<class STPDArray, class textOracle, class phiFunction>
class stpd_index{

private:

	phiFunction phi; // phi-function data structure
	textOracle O; // random access text oracle
	STPDArray S; // stpd array binary search

public:
	
	stpd_index(){} // empty constructor

	// standard index constructor
	void build_colex_m(const std::string &text_filepath, const std::string &sampling_filepath,
		               const std::string &rbwt_filepath, const std::string &pa_filepath, size_t refLen)
	{
		std::cout << "Constructing the STPD-index for " << text_filepath << std::endl;
		std::cout << "Step 1) Constructing the random-access text oracle..." << std::endl;
		if(refLen > 0){ O.build(text_filepath,refLen); }
		else{ O.build(text_filepath,1.0,0); }
		std::cout << "Step 2) Constructing the STPD-array binary search data structure..." << std::endl;
		S.build(sampling_filepath,&O,false); 
		std::cout << "Step 3) Constructing the phi function..." << std::endl;
	  	phi.build(rbwt_filepath,pa_filepath);
	  	
	  	std::cout << "Index succesfully built!" << std::endl;
	}
	// optimized index constructor
	void build_colex_m(const std::string &text_filepath, const std::string &sampling_filepath,
		               const std::string &rbwt_filepath, const std::string &pa_filepath,
		               const std::string &lcs_filepath, size_t refLen)
	{
		std::cout << "Constructing the STPD-index for " << text_filepath << std::endl;
		std::cout << "Step 1) Constructing the random-access text oracle..." << std::endl;
		if(refLen > 0){ O.build(text_filepath,refLen); }
		else{ O.build(text_filepath,1.0,0); }
		std::cout << "Step 2) Constructing the STPD-array binary search data structure..." << std::endl;
		S.build(text_filepath,sampling_filepath,lcs_filepath,pa_filepath,&O,false); 
		std::cout << "Step 3) Constructing the phi-function..." << std::endl;
	  	phi.build(rbwt_filepath,pa_filepath);
	  	
	  	std::cout << "Index succesfully built!" << std::endl;
	}

	/*
	void build_colex_pm(const std::string &text_filepath, const std::string &sampling_filepath,
		                const std::string &rbwt_filepath, const std::string &pa_filepath, size_t refLen)
	{
		std::cout << "Constructing the STPD-index for " << text_filepath << std::endl;
		std::cout << "Step 1) Constructing the random-access text oracle..." << std::endl;
		if(refLen > 0){ O.build(text_filepath,refLen); }
		else{ O.build(text_filepath,1.0,0); }
		std::cout << "Step 2) Constructing the STPD-array binary search data structure..." << std::endl;
		S.build(sampling_filepath,&O,true); 
		std::cout << "Step 3) Constructing the phi function..." << std::endl;
	  	phi.build(rbwt_filepath,pa_filepath);

	  	std::cout << "Index succesfully built!" << std::endl;
	}
	*/
	
	usafe_t store(const std::string &index_filepath)
	{
		std::ofstream out(index_filepath);

		usafe_t O_bytes = O.serialize(out);
		std::cout << "Storing random-access data structure, size = " << O_bytes << " bytes" << std::endl;
		usafe_t S_bytes = S.serialize(out);
		std::cout << "Storing STPD-array data structure, size = " << S_bytes << " bytes" << std::endl;
		usafe_t phi_size = phi.serialize(out);
		std::cout << "Storing phi data structure, size = " << phi_size << " bytes" << std::endl;

		std::cout << "## Index succesfully stored!" << std::endl;
		std::cout << "Total index size = " << O_bytes + S_bytes + phi_size << " bytes" << std::endl;
		
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

		std::cout << "Index succesfully loaded!" << std::endl;

		in.close();
	}

	// locate all occurrences exponential search
	std::tuple<std::vector<uint_t>,double,double> 
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
				return std::make_tuple(std::vector<uint_t>{},0,0);

			i_occ.second = std::get<0>(j);
			usafe_t f = O.LCP(pattern,i_occ.first,i_occ.second+1);
			i_occ.first = i_occ.first + f + 1;
			i_occ.second = i_occ.second + f;
		}
		std::chrono::duration<double> duration_mid = 
				std::chrono::high_resolution_clock::now() - start;

		usafe_t high = 2, low = 0;
		std::vector<uint_t> res{uint_t(i_occ.second)};
		while(true)
		{
			usafe_t phi_steps = high/2;
			while(phi_steps-- > 0)
			{
				i_occ.second = phi.phi_safe(i_occ.second);
				if(i_occ.second == -1)
				{
					high -= phi_steps;
					binary_search_occs(low,high,m,pattern,res);
					res.resize(low);

					std::chrono::duration<double> duration = 
							std::chrono::high_resolution_clock::now() - start;

					return std::make_tuple(res,duration.count(),duration_mid.count());			
				}

				res.push_back(i_occ.second);
			}

			usafe_t f = O.LCS(pattern,m-1,res[high-1]);
			if(f < m)
				break;

			low = high;
			high *= 2;
		}

		binary_search_occs(low,high,m,pattern,res);
		res.resize(low);

		std::chrono::duration<double> duration = 
				std::chrono::high_resolution_clock::now() - start;

		return std::make_tuple(res,duration.count(),duration_mid.count());
	}
	/*
	std::tuple<std::vector<uint_t>,double,double> 
						 locate_pattern(const std::string pattern) const
	{
		auto start = std::chrono::high_resolution_clock::now();

		usafe_t i = 1, m = pattern.size();
		int_t lower_occ, upper_occ;
		bool_t mismatch_found;

		this->lower_sample(pattern, lower_occ, mismatch_found);
		this->upper_sample(pattern, upper_occ);

		if(mismatch_found)
			return std::make_tuple(std::vector<uint_t>{},0,0);

		std::chrono::duration<double> duration_mid = 
				std::chrono::high_resolution_clock::now() - start;

		std::vector<uint_t> res{uint_t(lower_occ)};
		while(lower_occ != upper_occ)
		{
			lower_occ = phi.phi_unsafe(lower_occ);
			res.push_back(lower_occ);
		}
		
		std::chrono::duration<double> duration = 
				std::chrono::high_resolution_clock::now() - start;

		return std::make_tuple(res,duration.count(),duration_mid.count());
	}
	*/

	// run locate all occurrence queries on all patterns in a fasta file
	/*
		Parameters:
		- patternFile: FASTA file path containing the patterns
		- thr: Maximum number of occurrences to report for each pattern
		Output: A patternFile.occs file containing the positions of the 
		        patterns in the original text, and some statistics printed 
		        to the standard output

		Note that the check_occs_correctness function assumes that each pattern 
		occurs at least once in the text.
	*/
	void locate_fasta(const std::string patternFile, usafe_t thr) const
	{
		std::ifstream patterns(patternFile);
		std::ofstream   output(patternFile+".occs");

		std::string line, header;
		usafe_t i=0, c=0;
		std::tuple<std::vector<uint_t>,double,double> o;
		double tot_duration = 0;

		malloc_count_reset_peak();

		uint_t tot_occs = 0;
		while(std::getline(patterns, line))
		{
			if(i%2 != 0)
			{
				//if(this->S.is_index_large())
				//	o = locate_pattern(line);
				//else
					o = locate_pattern_exp_search(line);

				output << header << std::endl;
				if(std::get<0>(o).size() >= 0)
				{	
					usafe_t j=0;
					for(auto& e:std::get<0>(o))
					{
						output << e << " ";
						if(++j > thr-1) break;
					}
					output << std::endl;
				}

				tot_duration += std::get<1>(o);
				c += line.size();
				tot_occs += std::get<0>(o).size();
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
				  << "Total number of occurrences found = " << tot_occs << std::endl
		          << "Elapsed time per pattern = " <<
				     (tot_duration/(i/2))*1000000000 << " nanoSec" << std::endl
		          << "Elapsed time per character = " <<
				     (tot_duration/(c))*1000000000 << " nanoSec" << std::endl
		          << "Elapsed time per occurrence = " <<
				     (tot_duration/(tot_occs))*1000000000 << " nanoSec" << std::endl;
	}

	// check running time and correctness of locating all occurrences queries
	/*
		Parameters:
		- patternFile: FASTA file path containing the patterns
		Output: some statistics printed to the standard output
	*/
	void locate_fasta_test_running_time(const std::string patternFile) const
	{
		std::ifstream patterns(patternFile);

		std::string line, header;
		usafe_t i=0, c=0;
		std::tuple<std::vector<uint_t>,double,double> o;
		double tot_duration = 0, binary_search_duration = 0;

		malloc_count_reset_peak();

		uint_t tot_occs = 0;
		while(std::getline(patterns, line))
		{
			if(i%2 != 0)
			{
				//if(this->S.is_index_large())
				//	o = locate_pattern(line);
				//else
					o = locate_pattern_exp_search(line);

				tot_duration += std::get<1>(o);
				binary_search_duration += std::get<2>(o);
				c += line.size();
				tot_occs += std::get<0>(o).size();

				if(not check_occs_correctness(std::get<0>(o),line))
					exit(1);
			}
			else{ header = line; }
			i++;
		}

		patterns.close();

		std::cout << "Memory peak while running pattern matching queries = " <<
				     malloc_count_peak() << " bytes" << std::endl
		          << "Elapsed time while running pattern matching queries = " <<
				     tot_duration << " sec" << std::endl 
				  << "Total number of occurrences found = " << tot_occs << std::endl
		          << "Number of patterns = " << i/2 
		 		  << ", Total number of characters = " << c << std::endl
		          << "Elapsed time per pattern = " <<
				     (tot_duration/(i/2))*1000000000 << " nanoSec" << std::endl
		          << "Elapsed time per character = " <<
				     (tot_duration/(c))*1000000000 << " nanoSec" << std::endl
		          << "Elapsed time per occurrence = " <<
				     (tot_duration/(tot_occs))*1000000000 << " nanoSec" << std::endl
				  << "Elapsed time running binary search on the STPD samples vector = " <<
				  	 binary_search_duration << " sec" << std::endl
				  << "Elapsed time running the phi queries = " <<
				  	 tot_duration-binary_search_duration << " sec" << std::endl
				  << "Percentage time taken for running the phi queries = " <<
				     ((tot_duration-binary_search_duration)/tot_duration)*100 << "%" << std::endl;
	}

private:

	bool check_occs_correctness(const std::vector<uint_t>& occs, const std::string& patt) const
	{	
		if(occs.size() == 0)
		{
			std::cerr << "Pattern " << patt << " has no occurrences!" << std::endl;
			return false;
		}
		for(const auto& e : occs)
		{
			uint_t lcs = O.LCS(patt,patt.size()-1,e);
			if(lcs != patt.size())
			{
				std::cerr << "Error detected with pattern: " << patt << " and occ " << e << std::endl;
				return false;
			}
		}

		return true;
	}

	inline void binary_search_occs(usafe_t& low, usafe_t& high, usafe_t m, 
		                      const std::string& pattern, const std::vector<uint_t>& res) const
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
	
}; // stpd_index
}  // stpd

#endif