// Copyright (c) 2025, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 *  stpd-index: implementation of the Suffix Tree Path Decomposition index
 */

#pragma once

#include <chrono>
#include <malloc_count.h> 

#include <r-index_phi_inv_intlv.hpp>   
#include <move-r_phi_inv.hpp>             
#include <RLZ/RLZ_DNA_sux.hpp>            
#include <bitpacked_text_oracle.hpp>
#include <tabulated_binary_search_DNA.hpp>   

namespace stpd{

template<class binarySearch, class textOracle, class phiFunction>
class colex_minus_index{

private:

	phiFunction  P; // phi-function data structure
	textOracle   O;   // random access text oracle
	binarySearch S;   // sampled prefix array binary search

public:
	
	colex_minus_index(){} // empty constructor

	// optimized index constructor
	void build( const std::string &text_filepath, const std::string &sampling_filepath,
		        const std::string &rbwt_filepath, const std::string &pa_filepath,
		        const std::string &lcs_filepath,
		        const usafe_t reference_length, const usafe_t tabulation_length = 15 )
	{
		std::cout << "[INFO] Constructing the STPD-index using the path decomposition in " 
		          << sampling_filepath << "\n" << std::endl;
		std::cout << "[STEP 1] Constructing the random-access text oracle..." << std::endl;
		O.build(text_filepath,reference_length);
		std::cout << "[STEP 2] Constructing the binary search data structure..." << std::endl;
		S.build(text_filepath,sampling_filepath,lcs_filepath,pa_filepath,&O,tabulation_length); 
		std::cout << "[STEP 3] Constructing the phi function..." << "\n" << std::endl;
	  	P.build(rbwt_filepath,pa_filepath);
	  	
	  	std::cout << "[DONE] Index successfully built!" << "\n" << std::endl;
	}
	
	usafe_t store(const std::string &index_filepath)
	{
		std::ofstream out(index_filepath);
		std::cout << "[INFO] Writing components to disk:" << std::endl;

		usafe_t O_bytes  = O.serialize(out);
		std::cout << "		- Random-access data structure size = " << O_bytes << " bytes" << std::endl;
		usafe_t S_bytes  = S.serialize(out);
		std::cout << "		- Binary search data structure size = " << S_bytes << " bytes" << std::endl;
		usafe_t phi_size = P.serialize(out);
		std::cout << "		- Phi-function data structure size = " << phi_size << " bytes" << "\n" << std::endl;

		std::cout << "[DONE] Index successfully stored!" << std::endl;
		std::cout << "		→ Total index size in disk = " << O_bytes + S_bytes + phi_size << " bytes" << "\n" << std::endl;
		
		out.close();

		return O_bytes + S_bytes + phi_size;
	}
	
	void load(const std::string &index_filepath)
	{
		std::ifstream in(index_filepath);
		std::cout << "[INFO] Loading components to disk:" << std::endl;

		std::cout << "		- Random-access text oracle..." << std::endl;
		O.load(in);
		std::cout << "		- Binary search data structure..." << std::endl;
		S.load(in,&(this->O));
		std::cout << "		- Phi-function data structure..." << "\n" << std::endl;
		P.load(in);

		std::cout << "[DONE] Index successfully loaded!" << "\n" << std::endl;

		in.close();
	}

	// locate all occurrences exponential search
	std::vector<uint_t>
	locate_all_occurrences(const std::string &pattern,
						   const usafe_t     t = std::numeric_limits<usafe_t>::max(),
		                   const double      c = 1.5,
		                   const usafe_t     h = 2) const
	{
		usafe_t m = pattern.size();
		auto i_occ = S.match_first_prefix_with_tabulation(pattern);
		//std::cout << "first match -> " << i_occ.first << "," << i_occ.second << std::endl;
		//exit(1);

		if(i_occ.first-1 < m)
		{
			safe_t a, f;

			while(i_occ.first-1 < m)
			{
				a = S.binary_search_with_tabluation(pattern,0,i_occ.first);

				if( a < 0 ) 
					return std::vector<uint_t>{};

				f            = O.LCP(pattern,i_occ.first,a + 1);
				i_occ.first  = i_occ.first + f + 1;
				i_occ.second = a + f;
			}

			if( S.is_bad_anchor(pattern,i_occ.first-f-1,a) )
				return std::vector<uint_t>{};
		}

		usafe_t high = 1, low = 0;
		usafe_t j = t - 1;
		std::vector<uint_t> res{uint_t(i_occ.second)};
		P.init_phi(i_occ.second);

		while(j-- > 0)
		{
			safe_t phi_steps = high - low;
			while(phi_steps-- > 0)
			{
				i_occ.second = P.phi_next();
				res.push_back(i_occ.second);
			}

			usafe_t f = O.LCS(pattern,m-1,res[high-1]);

			if(f < m)
				break;

			low = high;
			high = std::ceil(high*c);
		}
		binary_search_occs(low,high,m,pattern,res);
		res.resize(low);

		return res;
	}
	
	// locate one occurrence 
	safe_t
	locate_first_occurrence(const std::string &pattern) const
	{
		usafe_t m  = pattern.size();
		auto i_occ = S.match_first_prefix_with_tabulation(pattern);

		if(i_occ.first-1 < m)
		{
			safe_t a, f;

			while(i_occ.first-1 < m)
			{
				a = S.binary_search_with_tabluation(pattern,0,i_occ.first);

				if( a < 0 ) 
					return -1;

				f            = O.LCP(pattern,i_occ.first,a+1);
				i_occ.first  = i_occ.first + f + 1;
				i_occ.second = a + f;
			}

			if( S.is_bad_anchor(pattern,i_occ.first-f-1,a) )
				return -1;
		}

		return i_occ.second;
	}

	// run locate all occurrence queries on all patterns in a fasta file
	/*
		Parameters:
		- pattern_file:  FASTA file path containing the patterns
		- occ_threshold: Maximum number of occurrences to report for each pattern
		Output: A patternFile.occs file containing the positions of the 
		        patterns in the original text, and some statistics printed 
		        to the standard output

		Note that the check_occs_correctness function assumes that each pattern 
		occurs at least once in the text.
	*/
	void
	locate_all_fasta(const std::string pattern_file, usafe_t occ_threshold,
		             bool_t check_correctness = false) const
	{
		malloc_count_reset_peak();
		auto start = std::chrono::high_resolution_clock::now();

		std::ifstream   patterns(pattern_file);
		std::ofstream   output  (pattern_file+".occs");

		std::string line, header;
		usafe_t i = 0, tot_char = 0;
		double tot_duration = 0;

		uint_t tot_occs = 0;
		while(std::getline(patterns, line))
		{
			if(i%2 != 0)
			{
				auto o = locate_all_occurrences(line, occ_threshold);

				if( check_correctness and (not are_occs_correct(o,line)) )
					exit(1);

				output << header << std::endl;
				if(o.size() >= 0)
				{	
					usafe_t j=0;
					for(auto& e : o)
						output << e << " ";
					output << std::endl;
				}

				tot_char += line.size();
				tot_occs += o.size();
			}
			else{ header = line; }
			i++;
		}

		patterns.close();
		output.close();

		std::chrono::duration<double> duration = 
				std::chrono::high_resolution_clock::now() - start;

		std::cout << "Memory peak while running pattern matching queries = " <<
				     malloc_count_peak() << " bytes" << std::endl
		          << "Elapsed time while running pattern matching queries = " <<
				     duration.count() << " sec" << std::endl 
		          << "Number of patterns = " << i/2 
		 		  << ", Total number of characters = " << tot_char << std::endl
				  << "Total number of occurrences found = " << tot_occs << std::endl;
	}

	void locate_one_fasta(const std::string patternFile, 
								bool_t check_correctness = false) const
	{
		malloc_count_reset_peak();
		auto start = std::chrono::high_resolution_clock::now();

		std::ifstream patterns(patternFile);
		std::ofstream   output(patternFile+".occs");

		std::string line, header;
		usafe_t i=0, tot_char=0;

		uint_t located = 0;
		while(std::getline(patterns, line))
		{
			if(i%2 != 0)
			{
				auto o = locate_first_occurrence(line);

				if( check_correctness and o == -1 ){
					std::cerr << "Unable to locate pattern: " 
							  << line << std::endl;
					exit(1);
				}
				if( check_correctness and
					(not are_occs_correct(std::vector<uint_t>{o},line)) ) 
					exit(1);
					
				output << header << std::endl;
				if(o != -1)
				{
					located++;
					output << o << std::endl;
				}

				tot_char += line.size();
			}
			else{ header = line; }
			i++;
		}

		patterns.close();
		output.close();

		std::chrono::duration<double> duration = 
				std::chrono::high_resolution_clock::now() - start;

		std::cout << "Memory peak while running pattern matching queries = " <<
				     malloc_count_peak() << " bytes" << std::endl
		          << "Elapsed time while running pattern matching queries = " <<
				     duration.count() << " sec" << std::endl 
		          << "Number of patterns = " << i/2 
		 		  << ", Total number of characters = " << tot_char << std::endl
				  << "Total number of successfully located patterns = " << located << std::endl;
	}

	void
	locate_all_benchmark(const std::string pattern_file, 
		                       usafe_t occ_threshold, 
		                       bool_t check_correctness = false) const
	{
		malloc_count_reset_peak();

		std::ifstream   patterns(pattern_file);
		std::ofstream   output  (pattern_file+".occs");
		std::ofstream   logfile (pattern_file+".locate_log");

		std::string line, header;
		usafe_t i = 0, tot_char = 0;
		double tot_duration = 0;

		uint_t tot_occs = 0;
		while(std::getline(patterns, line))
		{
			if(i%2 != 0)
			{
				auto start = std::chrono::high_resolution_clock::now();
				auto o = locate_all_occurrences(line, occ_threshold);
				std::chrono::duration<double> duration = 
						std::chrono::high_resolution_clock::now() - start;

				if( check_correctness and (not are_occs_correct(o,line)) )
					exit(1);

				output << header << std::endl;
				if(o.size() >= 0)
				{	
					usafe_t j=0;
					for(auto& e : o)
						output << e << " ";
					output << std::endl;
				}

				tot_duration += duration.count();
				tot_char += line.size();
				tot_occs += o.size();
			}
			else{ header = line; }
			i++;
		}

		patterns.close();
		output.close();

		std::cout << "Memory peak (bytes) = " << malloc_count_peak() << "\n"
		       << "Elapsed time (sec) = " << tot_duration << "\n"
		       << "Tot. occurences = " << tot_occs << "\n"
		       << "No. patterns = " << i/2 << "\n"
		       << "No. characters = " << tot_char << "\n"
		       << "Time per pattern (ns) = " << (tot_duration/(i/2))*1000000000 << "\n"
		       << "Time per character (ns) = " << (tot_duration/(tot_char))*1000000000 << "\n"
		       << "Time per occurrence (ns) = " << (tot_duration/(tot_occs))*1000000000 << "\n";

		logfile << "Memory peak (bytes) = " << malloc_count_peak() << "\n"
		        << "Elapsed time (sec) = " << tot_duration << "\n"
		        << "Tot. occurences = " << tot_occs << "\n"
		        << "No. patterns = " << i/2 << "\n"
		        << "No. characters = " << tot_char << "\n"
		        << "Time per pattern (ns) = " << (tot_duration/(i/2))*1000000000 << "\n"
		        << "Time per character (ns) = " << (tot_duration/(tot_char))*1000000000 << "\n"
		        << "Time per occurrence (ns) = " << (tot_duration/(tot_occs))*1000000000 << "\n";

		logfile.close();
	} 

	void locate_one_benchmark(const std::string pattern_file,
									bool_t check_correctness = false) const
	{
		malloc_count_reset_peak();

		std::ifstream patterns(pattern_file);
		std::ofstream output  (pattern_file+".occs");
		std::ofstream logfile (pattern_file+".primary_log");

		std::string line, header;
		usafe_t i=0, tot_char=0;
		double tot_duration = 0;

		uint_t located = 0;
		while(std::getline(patterns, line))
		{
			if(i%2 != 0)
			{
				auto start = std::chrono::high_resolution_clock::now();
				auto o = locate_first_occurrence(line);
				std::chrono::duration<double> duration = 
						std::chrono::high_resolution_clock::now() - start;

				if( o == -1 ){
					std::cerr << "Unable to locate pattern: " 
							  << line << std::endl;
					exit(1);
				}
				if( not are_occs_correct(std::vector<uint_t>{o},line) ) 
					exit(1);

				output << header << std::endl;
				if(o != -1)
				{
					located++;
					output << o << std::endl;
				}

				tot_duration += duration.count();
				tot_char += line.size();
			}
			else{ header = line; }
			i++;
		}

		patterns.close();
		output.close();

		std::cout << "Memory peak (bytes) = " << malloc_count_peak() << "\n"
		       << "Elapsed time (sec) = " << tot_duration << "\n"
		       << "No. patterns = " << i/2 << "\n"
		       << "No. characters = " << tot_char << "\n"
		       << "Time per pattern (ns) = " << (tot_duration/(i/2))*1000000000 << "\n"
		       << "Time per character (ns) = " << (tot_duration/(tot_char))*1000000000 << "\n";

		logfile << "Memory peak (bytes) = " << malloc_count_peak() << "\n"
		        << "Elapsed time (sec) = " << tot_duration << "\n"
		        << "No. patterns = " << i/2 << "\n"
		        << "No. characters = " << tot_char << "\n"
		        << "Time per pattern (ns) = " << (tot_duration/(i/2))*1000000000 << "\n"
		        << "Time per character (ns) = " << (tot_duration/(tot_char))*1000000000 << "\n";

		logfile.close();
	}

private:

	bool are_occs_correct(const std::vector<uint_t>& occs, const std::string& patt) const
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
				std::cerr << "Error detected with pattern: " << patt 
				          << " and occ " << e << std::endl;
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
	
}; // stpd_index
}  // stpd