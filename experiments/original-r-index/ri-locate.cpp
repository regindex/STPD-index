// Copyright (c) 2017, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#include <iostream>

#include <malloc_count.h>

#include "internal/r_index.hpp"

#include "internal/utils.hpp"

using namespace ri;
using namespace std;

string check = string();//check occurrences on this text
bool hyb=false;
string ofile;

void help(){
	cout << "ri-locate: locate all occurrences of the input patterns." << endl << endl;

	cout << "Usage: ri-locate [options] <index> <patterns>" << endl;
	cout << "   -c <text>    check correctness of each pattern occurrence on this text file (must be the same indexed)" << endl;
	//cout << "   -h           use hybrid bitvectors instead of elias-fano in both RLBWT and predecessor structures. -h is required "<<endl;
	//cout << "                if the index was built with -h options enabled."<<endl;
	cout << "   -o <ofile>   write pattern occurrences to this file (ASCII)" << endl;
	cout << "   <index>      index file (with extension .ri)" << endl;
	cout << "   <patterns>   file in FASTA format containing the patterns." << endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	if(s.compare("-c")==0){

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -c option." << endl;
			help();
		}

		check = string(argv[ptr]);
		ptr++;

	}else if(s.compare("-o")==0){

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -o option." << endl;
			help();
		}

		ofile = string(argv[ptr]);
		ptr++;

	}/*else if(s.compare("-h")==0){

		hyb=true;

	}*/else{

		cout << "Error: unknown option " << s << endl;
		help();

	}

}


template<class idx_t>
void locate(std::ifstream& in, string patterns){

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;

    string text;
    bool c = false;

    ofstream out;

    if(ofile.compare(string())!=0){

    	out = ofstream(ofile);

    }

    if(check.compare(string()) != 0){

    	c = true;

		ifstream ifs1(check);
		stringstream ss;
		ss << ifs1.rdbuf();//read the file
		text = ss.str();

    }

    auto t1 = high_resolution_clock::now();

    idx_t idx;

	idx.load(in);

	auto t2 = high_resolution_clock::now();

	cout << "searching patterns in " << patterns << endl;
	ifstream ifs(patterns);

	std::ofstream output(patterns+".occRindex");
	std::string line, header;
	size_t ca = 0;
	size_t tot_occs = 0;
	size_t i = 0;
	double tot_duration = 0;

	malloc_count_reset_peak();
	/*
	while(std::getline(ifs, line))
	{
		if(i%2 != 0)
		{
			auto OCC = idx.count_and_get_occ_time(line);

			output << header << std::endl;
			if(get<0>(OCC).second >= get<0>(OCC).first){ output << get<1>(OCC) << " " << line.size() << std::endl; }
			else{ output << "-1 " << line.size() << std::endl; }

			tot_duration += get<2>(OCC);
			ca += line.size();
		}
		else{ header = line; }
		i++;
	}
	*/
	while(std::getline(ifs, line))
	{
		if(i%2 != 0)
		{
			auto OCC = idx.locate_all_and_time(line);

			output << header << std::endl;
			if(OCC.first.size() > 0)
			{ 
				//for(const auto& e:OCC.first)
				//	output << e << " ";
				//output << std::endl;
				output << OCC.first.size() << endl;
			}
			else{ output << "-1 " << line.size() << std::endl; }
			
			tot_duration += OCC.second;
			ca += line.size();
			tot_occs += OCC.first.size();
		}
		else{ header = line; }
		i++;
	}

	ifs.close();
	output.close();

	auto t3 = high_resolution_clock::now();

	std::cout << "Memory peak while running pattern matching queries = " <<
		     malloc_count_peak() << " bytes" << std::endl
              << "Elapsed time while running pattern matching queries = " <<
		         tot_duration << " sec" << std::endl
		      << "Total number of occurrences found = " << tot_occs << std::endl
         	  << "Number of patterns = " << i/2 
 		  	  << ", Total number of characters = " << ca << std::endl
          	  << "Elapsed time per pattern = " <<
		         (tot_duration/(i/2))*1000000000  << " nanoSec" << std::endl
          	  << "Elapsed time per character = "  <<
		         (tot_duration/(ca))*1000000000   << " nanoSec" << std::endl
          	  << "Elapsed time per occurrence = " <<
		         (tot_duration/(tot_occs))*1000000000   << " nanoSec" << std::endl;

}

template<class idx_t>
void locate_test_running_time(std::ifstream& in, string patterns)
{
    string text;
    bool c = false;

    idx_t idx;

	idx.load(in);

	cout << "searching patterns in " << patterns << endl;
	ifstream ifs(patterns);

	std::string line, header;
	size_t ca = 0;
	size_t tot_occs = 0;
	size_t i = 0;
	double tot_duration = 0, backward_search_duration = 0;

	malloc_count_reset_peak();

	while(std::getline(ifs, line))
	{
		if(i%2 != 0)
		{
			auto OCC = idx.locate_all_and_time(line);
			
			tot_duration += std::get<1>(OCC);
			backward_search_duration += std::get<2>(OCC);
			ca += line.size();
			tot_occs += std::get<0>(OCC).size();
		}
		else{ header = line; }
		i++;
	}

	ifs.close();

	std::cout << "Memory peak while running pattern matching queries = " <<
		     malloc_count_peak() << " bytes" << std::endl
              << "Elapsed time while running pattern matching queries = " <<
		         tot_duration << " sec" << std::endl
		      << "Total number of occurrences found = " << tot_occs << std::endl
         	  << "Number of patterns = " << i/2 
 		  	  << ", Total number of characters = " << ca << std::endl
          	  << "Elapsed time per pattern = " <<
		         (tot_duration/(i/2))*1000000000  << " nanoSec" << std::endl
          	  << "Elapsed time per character = "  <<
		         (tot_duration/(ca))*1000000000   << " nanoSec" << std::endl
          	  << "Elapsed time per occurrence = " <<
		         (tot_duration/(tot_occs))*1000000000   << " nanoSec" << std::endl
			  << "Elapsed time running backward search algorithm = " <<
			  	 backward_search_duration << " sec" << std::endl
			  << "Elapsed time running phi queries = " <<
			  	 (tot_duration-backward_search_duration) << " sec" << std::endl
			  << "Percentage time taken for running the phi queries = " <<
			     ((tot_duration-backward_search_duration)/tot_duration)*100 << "%" << std::endl;
}

int main(int argc, char** argv){

	if(argc < 3)
		help();

	int ptr = 1;

	while(ptr<argc-2)
		parse_args(argv, argc, ptr);

	string idx_file(argv[ptr]);
	string patt_file(argv[ptr+1]);

	std::ifstream in(idx_file);

	bool fast;

	//fast or small index?
	in.read((char*)&fast,sizeof(fast));

	cout << "Loading r-index" << endl;

	if(hyb){

		//locate<r_index<sparse_hyb_vector,rle_string_hyb> >(in, patt_file);

	}else{

		locate_test_running_time<r_index<> >(in, patt_file);

	}

	in.close();

}
