// Copyright (c) 2025, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 *  pa-index: implementation of the Prefix Array index
 */

#ifndef PT_INDEX_HPP_
#define PT_INDEX_HPP_

#include <chrono>
#include <malloc_count.h> 

#include <common.hpp> 
#include <RLZ/RLZ_DNA_sux.hpp> 
#include <bitpacked_text_oracle.hpp>  

#include "prefix_tree.hpp"

namespace stpd{

template<class textOracle>
class pt_index{

private:

	textOracle O; 			  // random access text oracle
	PrefixTree<textOracle> S; // prefix tree data structure

public:
	
	pt_index(){} // empty constructor

	// standard index constructor
	void build(const std::string &text_filepath, const std::string &pa_filepath,
	        									 const std::string &lcs_filepath,
	        									            size_t reference_length)
	{
		std::cout << "[INFO] Constructing the PT index using the prefix array in " << pa_filepath << "\n" << std::endl;
		std::cout << "[STEP 1] Constructing the random-access text oracle..." << std::endl;
		O.build(text_filepath,reference_length);
		std::cout << "[STEP 2] Constructing the Prefix tree data structure..." << std::endl;
		S.build(&O,pa_filepath,lcs_filepath,false); 
	  	
	  	std::cout << "[DONE] Index successfully built!" << "\n" << std::endl;
	}
	
	usafe_t store(const std::string &index_filepath)
	{
		std::ofstream out(index_filepath);
		std::cout << "[INFO] Writing components to disk:" << std::endl;

		usafe_t O_bytes = O.serialize(out);
		std::cout << "		- Random-access data structure size = " << O_bytes << " bytes" << std::endl;
		usafe_t S_bytes = S.serialize(out);
		std::cout << "		- Prefix tree data structure size = " << S_bytes << " bytes" << std::endl;

		std::cout << "[DONE] Index successfully stored!" << std::endl;
		std::cout << "		→ Total index size in disk = " << O_bytes + S_bytes << " bytes" << "\n" << std::endl;
		
		out.close();

		return O_bytes + S_bytes;
	}
	
	void load(const std::string &index_filepath)
	{
		std::ifstream in(index_filepath);
		std::cout << "[INFO] Loading components to disk:" << std::endl;

		std::cout << "		- Random-access text oracle..." << std::endl;
		O.load(in);
		std::cout << "		- Prefix tree data structure..." << std::endl;
		S.load(in,&(this->O));

		std::cout << "[DONE] Index successfully loaded!" << "\n" << std::endl;

		in.close();
	}
	
	std::vector<uint_t>
		locate_pattern(const std::string pattern, usafe_t thr) const
	{
		auto start = std::chrono::high_resolution_clock::now();

		usafe_t i = 1, m = pattern.size();
		int_t lower_occ, upper_occ;

		Node* l = S.locate_locus(pattern);

		Node* b = S.locate_smallest_leaf(l);
		Node* e = S.locate_largest_leaf(l);

		std::vector<uint_t> res = S.get_SA_range_thr(b,e,thr);

		std::chrono::duration<double> duration = 
				std::chrono::high_resolution_clock::now() - start;
		
		return res;
	}
	

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
	void locate_fasta(const std::string patternFile, usafe_t thr) 
	{
		malloc_count_reset_peak();
		auto start = std::chrono::high_resolution_clock::now();

		std::ifstream patterns(patternFile);
		std::ofstream   output(patternFile+".PToccs");

		std::string line, header;
		usafe_t i=0, tot_char=0;
		std::vector<uint_t> res;

		uint_t tot_occs = 0;
		while(std::getline(patterns, line))
		{
			if(i%2 != 0)
			{
				res = locate_pattern(line, thr);

				output << header << std::endl;
				if(res.size() >= 0)
				{	
					usafe_t j=0;
					for(auto& e: res)
					{
						output << e << " ";
						if(++j > thr-1) break;
					}
					output << std::endl;
				}

				tot_char += line.size();
				tot_occs += res.size();
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

	// check running time and correctness of locating all occurrences queries
	/*
		Parameters:
		- patternFile: FASTA file path containing the patterns
		Output: some statistics printed to the standard output
	*/
	void locate_fasta_benchmark(const std::string patternFile, usafe_t thr) 
	{
		malloc_count_reset_peak();

		std::ifstream patterns(patternFile);
		std::ofstream   output(patternFile+".PToccs");
		std::ofstream  logfile(patternFile+".PT_locate_log");

		std::string line, header;
		usafe_t i=0, tot_char=0;
		std::vector<uint_t> res;
		double tot_duration = 0;

		uint_t tot_occs = 0;
		while(std::getline(patterns, line))
		{
			if(i%2 != 0)
			{
				auto start = std::chrono::high_resolution_clock::now();
				res = locate_pattern(line, thr);
				std::chrono::duration<double> duration = 
						std::chrono::high_resolution_clock::now() - start;

				output << header << std::endl;
				if(res.size() >= 0)
				{	
					usafe_t j=0;
					for(auto& e: res)
					{
						output << e << " ";
						if(++j > thr-1) break;
					}
					output << std::endl;
				}

				tot_duration += duration.count();
				tot_char += line.size();
				tot_occs += res.size();
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


}; // pt_index
}  // stpd

#endif