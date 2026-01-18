// Copyright (c) 2017, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/* Modified by Davide Cenzato to implement the STPD-index benchmarks */

#include <iostream>

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
	cout << "   -q <query_type>   select the type of query benchmark (loc_all|loc_prim|count)" << endl;
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
void locate_primary_benchmark(std::ifstream& in, string patterns)
{
    idx_t idx;
    cout << "Loading the r-index" << endl;
	idx.load(in);

	cout << "searching patterns in " << patterns << endl;
	ifstream ifs(patterns);
	std::ofstream output(patterns+".ri_occs");
	std::ofstream logfile(patterns+".ri_primary_log");

	std::string line, header;
	size_t i = 0, tot_char = 0, tot_occs = 0;
	double tot_duration = 0;

	malloc_count_reset_peak();

	while(std::getline(ifs, line))
	{
		if(i%2 != 0)
		{
			auto start = std::chrono::high_resolution_clock::now();
			auto RES = idx.count_and_get_occ(line);
			std::chrono::duration<double> duration = 
					std::chrono::high_resolution_clock::now() - start;

			output << header << std::endl;
			if( RES.first.second >= RES.first.first )
			{
				output << RES.second << "\n";
				tot_occs += RES.first.second - RES.first.first + 1;
			}
			else { output << "-1\n"; }
			
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
void locate_all_benchmark(std::ifstream& in, string patterns, uint64_t thr)
{
    idx_t idx;
    cout << "Loading the r-index" << endl;
	idx.load(in);

	cout << "searching patterns in " << patterns << endl;
	ifstream ifs(patterns);
	std::ofstream output(patterns+".ri_occs");
	std::ofstream logfile(patterns+".ri_locate_log");

	std::string line, header;
	size_t i = 0, tot_char = 0, tot_occs = 0;
	double tot_duration = 0;

	malloc_count_reset_peak();

	while(std::getline(ifs, line))
	{
		if(i%2 != 0)
		{
			auto start = std::chrono::high_resolution_clock::now();
			auto OCC = idx.locate_all(line, thr);
			std::chrono::duration<double> duration = 
					std::chrono::high_resolution_clock::now() - start;

			output << header << std::endl;
			if(OCC.size() > 0)
			{ 
				for( const auto& e : OCC )
					output << e << " ";
				output << std::endl;
			}
			
			tot_duration += duration.count();
			tot_char += line.size();
			tot_occs += OCC.size();
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

	output.close();
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

	if( qtype == "loc_all" )
		locate_all_benchmark< r_index<> >(in, patt_file, maxOcc);
	else if( qtype == "loc_prim" )
		locate_primary_benchmark< r_index<> >(in, patt_file);
	else if( qtype == "count" )
		count_benchmark< r_index<> >(in, patt_file);
	else{
		cout << "Unspecified query type, see the help function (-t flag). Exiting..." << endl;
		exit(1);
	}

	in.close();
}