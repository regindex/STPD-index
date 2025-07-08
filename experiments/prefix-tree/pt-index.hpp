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

#include <RLZ_DNA_sux.hpp> // rlz random access text orcale
#include <common.hpp>.     // common definitions
#include "prefix_tree.hpp" // binary search ds

namespace stpd{

template<class textOracle>
class pt_index{

private:

	textOracle O; // random access text oracle
	PrefixTree<> S; // prefix tree data structure

public:
	
	pt_index(){} // empty constructor

	// standard index constructor
	void build(const std::string &text_filepath, const std::string &pa_filepath,
	        									 const std::string &lcs_filepath, size_t refLen)
	{
		std::cout << "[INFO] Constructing the PT index using the prefix array in " << pa_filepath << "\n" << std::endl;
		std::cout << "[STEP 1] Constructing the random-access text oracle..." << std::endl;
		if(refLen > 0){ O.build(text_filepath,refLen); }
		else{ O.build(text_filepath,1.0,0); }
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
	
	std::tuple<std::vector<uint_t>,double,double> 
						 locate_pattern(const std::string pattern) const
	{
		auto start = std::chrono::high_resolution_clock::now();

		usafe_t i = 1, m = pattern.size();
		int_t lower_occ, upper_occ;
		bool_t mismatch_found;

		Node* l = S.locate_locus(pattern);

		Node* b = S.locate_smallest_leaf(l);
		Node* e = S.locate_largest_leaf(l);

		std::chrono::duration<double> duration_mid = 
				std::chrono::high_resolution_clock::now() - start;

		std::vector<uint_t> res = S.get_SA_range(b,e);

		std::chrono::duration<double> duration = 
				std::chrono::high_resolution_clock::now() - start;
		
		return std::make_tuple(res,duration.count(),duration_mid.count());
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
		std::ifstream patterns(patternFile);
		std::ofstream   output(patternFile+".PAoccs");

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
				o = locate_pattern(line);

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
	void locate_fasta_test_running_time(const std::string patternFile) 
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
				//line = "CATCAACA";
				//std::cout << "pattern = " << line << std::endl;
				std::cout << "pattern = " << line << " -> ";

				o = locate_pattern(line);

				tot_duration += std::get<1>(o);
				binary_search_duration += std::get<2>(o);
				c += line.size();
				tot_occs += std::get<0>(o).size();

				std::cout << std::get<0>(o).size() << std::endl;

				if(not check_occs_correctness(std::get<0>(o),line))
					exit(1);

				//exit(1);
				//break;
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


}; // pt_index
}  // stpd

#endif