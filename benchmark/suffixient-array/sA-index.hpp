// Copyright (c) 2025, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 *  suffixient_index: implementation of the suffixient array-based index
 */

#pragma once

#include <chrono>
#include <malloc_count.h> 

#include <r-index_phi_inv.hpp>   
#include <move-r_phi_inv.hpp>    
#include <RLZ/RLZ_DNA_sux.hpp>            
#include <bitpacked_text_oracle.hpp>
#include <tabulated_binary_search_DNA.hpp>  

namespace suffixient{

	template<class suffArray, class textOracle>
	class suffixient_array_index{

	private:
		// text oracle data structure
		textOracle O;
		// suffixient array search data structure
		suffArray S;

	public:
		// empty constructor
		suffixient_array_index(){}

		// build suffixient index by indexing the supermaximal extensions
		void build(const std::string &text_filepath,
				   const std::string &sampling_filepath,
				   const     usafe_t reference_length,
				   const     usafe_t tabulation_length = 15)
		{
			std::cout << "[STEP 1] Constructing the random-access text oracle for "
			          << text_filepath << std::endl;
		  	O.build(text_filepath,reference_length);
			std::cout << "[STEP 2] Constructing the binary search data structure for "
					  << sampling_filepath << std::endl;
			S.build(text_filepath, sampling_filepath, &O, tabulation_length);
		}

		usafe_t store(const std::string &index_filepath)
		{
			std::ofstream out(index_filepath);

			usafe_t w_bytes = 0;
			w_bytes += O.serialize(out);
			w_bytes += S.serialize(out);
			
			out.close();

			return w_bytes;
		}

		void load(const std::string &index_filepath)
		{
			std::ifstream in(index_filepath);

			O.load(in);
			S.load(in,&(this->O));

			in.close();
		} 

		// locate primary occurrence 
		safe_t
		locate_primary_occurrence(const std::string &pattern) const
		{
			usafe_t m = pattern.size();
			auto i_occ = S.match_first_prefix_with_tabulation_suffixient(pattern);

			if(i_occ.first-1 < m)
			{
				safe_t a, f;

				while(i_occ.first-1 < m)
				{
					a = S.binary_search_with_tabluation_suffixient(pattern,0,i_occ.first);

					if( a < 0 ) 
						return -1;

					f            = O.LCP(pattern,i_occ.first,a + 1);
					i_occ.first  = i_occ.first + f + 1;
					i_occ.second = a + f;
				}

				if( S.is_bad_anchor(pattern,i_occ.first-f-1,a) )
					return -1;
			}

			return i_occ.second;
		}

		void locate_primary_benchmark(const std::string pattern_file) const
		{
		    // --- RESET MEMORY COUNTER ---
		    // compute index size + basic program overhead.
		    size_t index_memory_bytes = malloc_count_current();

		    // --- PHASE 1: LOAD DATA ---
		    // read patterns into memory first to avoid disk I/O latency.
		    std::ifstream patterns(pattern_file);
		    std::ofstream logfile(pattern_file + ".locate_log");

		    struct QueryData {
		        std::string header;
		        std::string pattern;
		    };
		    std::vector<QueryData> queries;
		    std::vector<safe_t> primary_occurrences;

		    std::string line, header;
		    usafe_t i = 0, tot_char = 0;

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

		    usafe_t tot_pattern = queries.size();
		    double  tot_primary_duration = 0;
		    usafe_t tot_primary_found = 0;

		    // --- PHASE 2: (Primary occurrences) ---
		    for(const auto& q : queries)
		    {
		        auto start = std::chrono::high_resolution_clock::now();
		        safe_t pocc = locate_primary_occurrence(q.pattern);
		        std::chrono::duration<double> duration = 
		                std::chrono::high_resolution_clock::now() - start;

		        tot_primary_duration += duration.count();
		        tot_primary_found += (pocc == -1 ? 0 : 1);
		    }

			// --- PHASE 4: LOG DATA ---
		    std::cout << "No. primary occs found = "    << tot_primary_found                    << "\n"
		              << "No. processed patterns = "    << tot_pattern                          << "\n"
		              << "No. processed characters = "  << tot_char                             << "\n" 
		              << "Index size (bytes) = "        << index_memory_bytes                   << "\n"
		              << "Time to find the primary occ. (sec) = "     << tot_primary_duration   << "\n"
		              << "Time to process one character (ns) = "      << 
		              						    (tot_primary_duration/tot_char) * 1000000000    << "\n"
		              << "Time per primary occurrence (ns) = "        << 
		                                        (tot_primary_duration/tot_pattern) * 1000000000 << std::endl;

		    logfile   << "No. primary occs found = "    << tot_primary_found                    << "\n"
		              << "No. processed patterns = "    << tot_pattern                          << "\n"
		              << "No. processed characters = "  << tot_char                             << "\n" 
		              << "Index size (bytes) = "        << index_memory_bytes                   << "\n"
		              << "Time to find the primary occ. (sec) = "     << tot_primary_duration   << "\n"
		              << "Time to process one character (ns) = "      << 
		              						    (tot_primary_duration/tot_char) * 1000000000    << "\n"
		              << "Time per primary occurrence (ns) = "        << 
		                                        (tot_primary_duration/tot_pattern) * 1000000000 << std::endl;

		    logfile.close();
		}
	};
}