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

	std::pair< std::pair<Node*,Node*>, safe_t >
	compute_range_and_occ(const std::string pattern)
	{
		Node* l = S.locate_locus(pattern);

		if( l == nullptr )
			return std::make_pair( std::make_pair(nullptr,nullptr), -1 );

		Node* b = S.locate_smallest_leaf(l);
		Node* e = S.locate_largest_leaf(l);

		return std::make_pair( std::make_pair(b,e), S[b] );
	}

	std::vector<uint_t>
	enumerate_secondary_occs(std::pair<Node*,Node*> range, usafe_t thr)
	{
		std::vector<uint_t> res {};
		
		Node* b = range.first;
		Node* e = range.second;

		if( b != e )
		{
			b = b->next_leaf;
			res = S.get_SA_range_thr(b,e,thr);
		}

		return res;
	}
	
	std::vector<uint_t>
		locate_pattern(const std::string pattern, usafe_t thr) const
	{
		Node* l = S.locate_locus(pattern);

		Node* b = S.locate_smallest_leaf(l);
		Node* e = S.locate_largest_leaf(l);

		std::vector<uint_t> res = S.get_SA_range_thr(b,e,thr);
		
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
	    // compute index size + basic program overhead.
	    size_t index_memory_bytes = malloc_count_current();

	    // --- PHASE 1: LOAD DATA ---
		// read patterns into memory first to avoid disk I/O latency.
		std::ifstream patterns(patternFile);
		std::ofstream   output(patternFile+".PToccs");
		std::ofstream  logfile(patternFile+".PT_locate_log");

	    struct QueryData {
	        std::string header;
	        std::string pattern;
	    };
	    std::vector<QueryData> queries;

	    using Step1ReturnType = decltype(compute_range_and_occ(""));
	    std::vector<Step1ReturnType> step1_results;

	    std::string line, header;
	    uint64_t i = 0, tot_char = 0;

	    while(std::getline(patterns, line)) 
	    {
	        if(i % 2 != 0) {
	            queries.push_back({header, line});
	            tot_char += line.size();
	        } else {
	            header = line;
	        }
	        i++;
	    }
	    patterns.close();

		// reserve space to prevent reallocation during timing
	    step1_results.reserve(queries.size());

	    uint64_t tot_pattern = queries.size();
	    uint64_t tot_soccs = 0;
	    double tot_primary_duration = 0;
	    double tot_secondary_duration = 0;

	    // --- PHASE 2: RUN STEP 1 (Primary occurrences) ---
	    for(const auto& q : queries)
	    {
	        auto start = std::chrono::high_resolution_clock::now();
			auto res = compute_range_and_occ(q.pattern);
	        std::chrono::duration<double> duration = 
	                std::chrono::high_resolution_clock::now() - start;

	        tot_primary_duration += duration.count();
	        step1_results.push_back(res);
	    }

	    // --- RESET MEMORY COUNTER ---
	    size_t memory_after_step_1 = malloc_count_current();
	    malloc_count_reset_peak();

	    // --- PHASE 3: RUN STEP 2 (Secondary occurrences) ---
	    // Step 2 enumerates the secondary occurences.
	    for(size_t k = 0; k < queries.size(); ++k)
	    {
	        auto& res = step1_results[k]; 
	        int64_t pocc = res.second;

	        auto start = std::chrono::high_resolution_clock::now();
	        auto soccs = enumerate_secondary_occs(res.first, thr - 1);
	        std::chrono::duration<double> duration = 
	                std::chrono::high_resolution_clock::now() - start;

	        tot_secondary_duration += duration.count();
	        tot_soccs += soccs.size();

	        // write output outside the timer
	        output << queries[k].header << std::endl;
	        if( pocc != -1 )
	        {   
	            output << pocc << " ";
	            for(auto& occ : soccs)
	                output << occ << " ";
	            output << std::endl;
	        }
	        else { output << std::endl; }
	    }

	    output.close();

		// --- PHASE 4: LOG DATA ---
	    std::cout << "No. occurences found = "      << tot_soccs + tot_pattern              << "\n"
	              << "No. processed patterns = "    << tot_pattern                          << "\n"
	              << "No. processed characters = "  << tot_char                             << "\n" 
	              << "Index size (bytes) = "        << index_memory_bytes                   << "\n"
	              << "Memory peak secondary occ.s (bytes) = "     << malloc_count_peak() -
	              							  (memory_after_step_1 - index_memory_bytes)    << "\n"
	              << "Time to find the primary occ. (sec) = "     << tot_primary_duration   << "\n"
	              << "Time to find the secondary occ.s (sec) = "  << tot_secondary_duration << "\n"
	              << "Time per primary occurrence (ns) = "        << 
	                                         (tot_primary_duration/tot_pattern) * 1000000000 << "\n"
	              << "Time per secondary occurrence (ns) = "      << 
	                                         (tot_secondary_duration/tot_soccs) * 1000000000 << std::endl;

	    logfile   << "No. occurences found = "      << tot_soccs + tot_pattern              << "\n"
	              << "No. processed patterns = "    << tot_pattern                          << "\n"
	              << "No. processed characters = "  << tot_char                             << "\n" 
	              << "Index size (bytes) = "        << index_memory_bytes                   << "\n"
	              << "Memory peak secondary occ.s (bytes) = "     << malloc_count_peak() -
	              							  (memory_after_step_1 - index_memory_bytes)    << "\n"
	              << "Time to find the primary occ. (sec) = "     << tot_primary_duration   << "\n"
	              << "Time to find the secondary occ.s (sec) = "  << tot_secondary_duration << "\n"
	              << "Time per primary occurrence (ns) = "        << 
	                                         (tot_primary_duration/tot_pattern) * 1000000000 << "\n"
	              << "Time per secondary occurrence (ns) = "      << 
	                                         (tot_secondary_duration/tot_soccs) * 1000000000 << std::endl;

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