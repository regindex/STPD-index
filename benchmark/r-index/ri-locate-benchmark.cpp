// Copyright (c) 2017, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/* Modified by Davide Cenzato to implement the STPD-index benchmarks */

#include <iostream>
#include <iomanip>

#include <malloc_count.h>

#include "internal/r_index.hpp"

#include "internal/utils.hpp"

using namespace ri;
using namespace std;

string qtype;
uint64_t maxOcc = (1ULL << 63) | ((1ULL << 63) - 1);

void help(){
	cout << "ri-locate: run locate/count query benchmarks." << endl << endl;

	cout << "Usage: ri-locate [options] <index> <patterns>" << endl;
	cout << "   -q <query_type>   select the type of query benchmark (secondary|primary)" << endl;
    cout << "   -t <max_occs>     maximum number of occurrences to report per pattern. (Def. all)" << endl;
	cout << "   <index>           index file (with extension .ri)" << endl;
	cout << "   <patterns>        file in FASTA format containing the patterns." << endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	if( s.compare("-t") == 0 ) {

		if(ptr>=argc-1) {
			cout << "Error: missing parameter after -t option." << endl;
			help();
		}

		maxOcc = stoull(argv[ptr]);
		ptr++;
	}
	else if( s.compare("-q") == 0 ){

		if( ptr >= argc-1 ) {
			cout << "Error: missing parameter after -q option." << endl;
			help();
		}

		qtype = string(argv[ptr]);
		ptr++;
	}
	else{

		cout << "Error: unknown option " << s << endl;
		help();
	}
}

template<class idx_t>
void count_benchmark(std::ifstream& in, string patterns)
{
    idx_t idx;
    cout << "Loading the r-index" << endl;
	idx.load(in);

	cout << "counting patterns in " << patterns << endl;
	ifstream ifs(patterns);
	std::ofstream output(patterns+".noccs");
	std::ofstream logfile(patterns+".count_log");

	std::string line, header;
	size_t i = 0, tot_char = 0, tot_occs = 0;
	double tot_duration = 0;

	malloc_count_reset_peak();

	while(std::getline(ifs, line))
	{
		if(i%2 != 0)
		{
			auto start = std::chrono::high_resolution_clock::now();
			auto range = idx.count(line);
			std::chrono::duration<double> duration = 
					std::chrono::high_resolution_clock::now() - start;

			output << header << std::endl;
			if( range.second >= range.first )
			{ 
				output << range.second - range.first - 1 << "\n";
				tot_occs += range.second - range.first - 1;
			}
			else { output << "0\n"; }
			
			tot_duration += duration.count();
			tot_char += line.size();
		}
		else{ header = line; }
		i++;
	}

	ifs.close();

	std::cout << "Memory peak (bytes) = " << malloc_count_peak() << "\n"
	       	  << "Elapsed time (sec) = " << tot_duration << "\n"
	       	  << "Tot. occurences = " << tot_occs << "\n"
	       	  << "No. patterns = " << i/2 << "\n"
	       	  << "No. characters = " << tot_char << "\n"
	       	  << "Time per pattern (ns) = " << (tot_duration/(i/2))*1000000000 << "\n"
	       	  << "Time per character (ns) = " << (tot_duration/(tot_char))*1000000000 << "\n";

	logfile << "Memory peak (bytes) = " << malloc_count_peak() << "\n"
	        << "Elapsed time (sec) = " << tot_duration << "\n"
	        << "Tot. occurences = " << tot_occs << "\n"
	        << "No. patterns = " << i/2 << "\n"
	        << "No. characters = " << tot_char << "\n"
	        << "Time per pattern (ns) = " << (tot_duration/(i/2))*1000000000 << "\n"
	        << "Time per character (ns) = " << (tot_duration/(tot_char))*1000000000 << "\n";

	output.close();
	logfile.close();
}

template<class idx_t>
void locate_benchmark(std::ifstream& in, string patterns, uint64_t occ_thr)
{
    idx_t idx;
    cout << "Loading the r-index" << endl;
	idx.load(in);

    // --- RESET MEMORY COUNTER ---
    // compute index size + basic program overhead.
    size_t index_memory_bytes = malloc_count_current();

    // --- PHASE 1: LOAD DATA ---
	// read patterns into memory first to avoid disk I/O latency.
	std::ifstream ifs(patterns);
	std::ofstream output(patterns+".ri_occs");
	std::ofstream logfile(patterns+".ri_primary_log");

    struct QueryData {
        std::string header;
        std::string pattern;
    };
    std::vector<QueryData> queries;

    using Step1ReturnType = decltype(idx.count_and_get_occ(""));
    std::vector<Step1ReturnType> step1_results;

    std::string line, header;
    uint64_t i = 0, tot_char = 0;

    while(std::getline(ifs, line)) 
    {
        if(i % 2 != 0) {
            queries.push_back({header, line});
            tot_char += line.size();
        } else {
            header = line;
        }
        i++;
    }
    ifs.close();

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
		auto res = idx.count_and_get_occ(q.pattern);
		res.first.first += 1;
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
        auto soccs = idx.locate_occurrences_in_range(res.first, pocc,
        											 occ_thr - 1);
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

template<class idx_t>
void locate_primary_benchmark(std::ifstream& in, string patterns)
{
    idx_t idx;
    cout << "Loading the r-index" << endl;
	idx.load(in,true);

    // --- RESET MEMORY COUNTER ---
    // compute index size + basic program overhead.
    size_t index_memory_bytes = malloc_count_current();

    // --- PHASE 1: LOAD DATA ---
	// read patterns into memory first to avoid disk I/O latency.
	std::ifstream ifs(patterns);
	std::ofstream logfile(patterns+".ri_primary_log");

    struct QueryData {
        std::string header;
        std::string pattern;
    };
    std::vector<QueryData> queries;

    std::string line, header;
    uint64_t i = 0, tot_char = 0;

    while(std::getline(ifs, line)) 
    {
        if(i % 2 != 0) {
            queries.push_back({header, line});
            tot_char += line.size();
        } else {
            header = line;
        }
        i++;
    }
    ifs.close();

    uint64_t tot_pattern = queries.size();
    double  tot_primary_duration = 0;
    uint64_t tot_primary_found = 0;

    // --- PHASE 2: (Primary occurrences) ---
    for(const auto& q : queries)
    {
        auto start = std::chrono::high_resolution_clock::now();
		auto res = idx.count_and_get_occ(q.pattern);
		int64_t pocc = res.second;
        std::chrono::duration<double> duration = 
                std::chrono::high_resolution_clock::now() - start;

        tot_primary_duration += duration.count();
        tot_primary_found += (res.first.first > res.first.second ? 0 : 1);
    }

	// --- PHASE 4: LOG DATA ---
    std::cout << "No. primary occs found = "    << tot_primary_found                    << "\n"
              << "No. processed patterns = "    << tot_pattern                          << "\n"
              << "No. processed characters = "  << tot_char                             << "\n" 
              << "Index size (bytes) = "        << index_memory_bytes                   << "\n"
              << "Time to find the primary occ. (sec) = "     << tot_primary_duration   << "\n"
              << "Time per primary occurrence (ns) = " << std::fixed << std::setprecision(2) << 
                                        (tot_primary_duration/tot_pattern) * 1000000000 << std::endl;

    logfile   << "No. primary occs found = "    << tot_primary_found                    << "\n"
              << "No. processed patterns = "    << tot_pattern                          << "\n"
              << "No. processed characters = "  << tot_char                             << "\n" 
              << "Index size (bytes) = "        << index_memory_bytes                   << "\n"
              << "Time to find the primary occ. (sec) = "     << tot_primary_duration   << "\n"
              << "Time per primary occurrence (ns) = " << std::fixed << std::setprecision(2) << 
                                        (tot_primary_duration/tot_pattern) * 1000000000 << std::endl;

    logfile.close();
}

template<class idx_t>
void locate_secondary_benchmark(std::ifstream& in, string patterns, uint64_t occ_thr)
{
    idx_t idx;
    cout << "Loading the r-index" << endl;
	idx.load(in);

    // --- RESET MEMORY COUNTER ---
    // compute index size + basic program overhead.
    size_t index_memory_bytes = malloc_count_current();

    // --- PHASE 1: LOAD DATA ---
	// read patterns into memory first to avoid disk I/O latency.
	std::ifstream ifs(patterns);
	std::ofstream logfile(patterns+".ri_primary_log");

    struct QueryData {
        std::string header;
        std::string pattern;
    };
    std::vector<QueryData> queries;

    std::string line, header;
    uint64_t i = 0, tot_char = 0;

    while(std::getline(ifs, line)) 
    {
        if(i % 2 != 0) {
            queries.push_back({header, line});
            tot_char += line.size();
        } else {
            header = line;
        }
        i++;
    }
    ifs.close();

    malloc_count_reset_peak();

    uint64_t tot_pattern = queries.size();
    uint64_t tot_soccs = 0;
    double tot_secondary_duration = 0;

    // --- PHASE 2: RUN QUERIES (Secondary occurrences) ---
    for(const auto& q : queries)
    {
        auto start = std::chrono::high_resolution_clock::now();
		auto res = idx.count_and_get_occ(q.pattern);
		res.first.first += 1;

		std::vector<uint32_t> soccs{};
		if(res.second != -1)
			soccs = idx.locate_occurrences_in_range(res.first, res.second,
        											occ_thr - 1);

        std::chrono::duration<double> duration = 
                std::chrono::high_resolution_clock::now() - start;

        tot_secondary_duration += duration.count();
        tot_soccs += soccs.size() + (res.second == -1 ? 0 : 1);
    }

	// --- PHASE 4: LOG DATA ---
    std::cout << "No. occurences found = "      << tot_soccs                            << "\n"
              << "No. processed patterns = "    << tot_pattern                          << "\n"
              << "No. processed characters = "  << tot_char                             << "\n" 
              << "Index size (bytes) = "        << index_memory_bytes                   << "\n"
              << "Memory peak secondary occ.s (bytes) = "     << malloc_count_peak() 	<< "\n"
              << "Time to find the secondary occ.s (sec) = "  << tot_secondary_duration << "\n"
              << "Time per secondary occurrence (ns) = "      << std::fixed << std::setprecision(2) <<  
                                         (tot_secondary_duration/tot_soccs) * 1000000000 << std::endl;

    logfile   << "No. occurences found = "      << tot_soccs + tot_pattern              << "\n"
              << "No. processed patterns = "    << tot_pattern                          << "\n"
              << "No. processed characters = "  << tot_char                             << "\n" 
              << "Index size (bytes) = "        << index_memory_bytes                   << "\n"
              << "Memory peak secondary occ.s (bytes) = "     << malloc_count_peak() 	<< "\n"
              << "Time to find the secondary occ.s (sec) = "  << tot_secondary_duration << "\n"
              << "Time per secondary occurrence (ns) = "      << std::fixed << std::setprecision(2) << 
                                         (tot_secondary_duration/tot_soccs) * 1000000000 << std::endl;

    logfile.close();
}

int main(int argc, char** argv){

	if(argc < 3)
		help();

	int ptr = 1;

	while( ptr < argc-2 )
		parse_args(argv, argc, ptr);

	string idx_file(argv[ptr]);
	string patt_file(argv[ptr+1]);

	std::ifstream in(idx_file);

	bool fast;
	//fast or small index?
	in.read((char*)&fast,sizeof(fast));

	if( qtype == "secondary" )
	{
		uint64_t maxOcc = (1ULL << 63) | ((1ULL << 63) - 1);
		locate_secondary_benchmark< r_index<> >(in, patt_file, maxOcc);
	}
	else if( qtype == "primary" )
		locate_primary_benchmark< r_index<> >(in, patt_file);
	else{
		cout << "Unspecified query type, see the help function (-t flag). Exiting..." << endl;
		exit(1);
	}

	in.close();
}