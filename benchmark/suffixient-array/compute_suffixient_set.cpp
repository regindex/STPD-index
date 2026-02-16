// Copyright (c) 2024, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#include <iostream>
#include <sdsl/construct.hpp>
#include <set>
#include <limits>
#include <algorithm>

using namespace std;
using namespace sdsl;

int_vector<8> T;
int_vector_buffer<> SA;
int_vector_buffer<> LCP;

size_t number_of_insertions = 0;

inline uint8_t BWT(uint64_t i){
	return SA[i] == 0 ? 0 : T[SA[i] - 1];
}

struct lcp_maxima{
	int64_t len;
	uint64_t pos;
	bool active;
};

void help(){

	cout << "suffixient [options]" << endl <<
	"Input: non-empty ASCII file without character 0x0, from file. Output: smallest suffixient set." << endl <<
	"Warning: if 0x0 appears, the standard input is read only until the first occurrence of 0x0 (excluded)." << endl <<
	"Options:" << endl <<
	"-h          Print usage info." << endl << 
	"-i <arg>    Input text file." << endl <<
	"-o <arg>    Store output to file using 40-bits unsigned integers. If not specified, output is streamed to standard output in human-readable format." << endl <<
	"-s          Sort output. Default: false." << endl <<
	"-p          Print to standard output size of suffixient set. Default: false." << endl <<
	"-r          Print to standard output number of equal-letter runs in the BWT of reverse text. Default: false." << endl;
	exit(0);
} 

inline void eval(uint64_t sigma,
								  int64_t m, 
								 vector<lcp_maxima>& R,
								 vector<vector<uint64_t>>& S)
{
	for(uint8_t c = 1; c < sigma; ++c)
	  if(m < R[c].len)
		  {
		    // process an active candidate
		    if(R[c].active)
		    {
		    	number_of_insertions += 1;
		    	S[c].push_back(R[c].pos - 1);
		    }

		    // update to inactive state
		    R[c] = {m,0,false};
		  }
}

int main(int argc, char** argv) {

	if(argc < 2) help();

	string input_file, output_file;

	bool sort = false;
	bool rho = false;
	bool runs = false;

	int opt;
	while ((opt = getopt(argc, argv, "prsho:i:")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			case 'o':
				output_file = string(optarg);
			break;
			case 'i':
				input_file = string(optarg);
			break;
			case 's':
				sort=true;
			break;
			case 'p':
				rho=true;
			break;
			case 'r':
				runs=true;
			break;
			default:
				help();
			return -1;
		}
	}

	cache_config cc;
	uint64_t N = 0; //including 0x0 terminator
	uint8_t sigma = 128; // alphabet size (including terminator 0x0)
	int64_t m = std::numeric_limits<int64_t>::max();
	uint64_t bwtruns = 1;

	{
		int_vector< 8> T_;
		load_vector_from_file( T_, input_file, 1 );
		append_zero_symbol( T_ );
		N = T_.size();

		if(N<2){
			cerr << "Error: empty text" <<  endl;
			help();
		}

		T = int_vector<8>(N - 1);

		for(uint64_t i = 0; i < N - 1; ++i)
		{
			uint8_t c = T_[N - i - 2];
			T[i] = c;
		}

		append_zero_symbol(T);
		store_to_cache(T, conf::KEY_TEXT, cc);
		construct_sa<8>(cc);
		construct_lcp_kasai<8>(cc);
		SA = int_vector_buffer<>(cache_file_name(conf::KEY_SA, cc));
		LCP = int_vector_buffer<>(cache_file_name(conf::KEY_LCP, cc));
	}

	vector<lcp_maxima> R(sigma,{-1,0,false}); //vector with candidate suffixient right-extensions
	vector<vector<uint64_t>> S(sigma,vector<uint64_t>{});

	for(uint64_t i=1;i<N;++i)
	{
		m = std::min(m,int64_t(LCP[i]));

		if(BWT(i) != BWT(i-1))
		{
			eval(sigma,m,R,S);

			for(uint64_t ip = i-1; ip < i+1; ++ip)
				if(int64_t(LCP[i]) > R[BWT(ip)].len)
					R[BWT(ip)] = {int64_t(LCP[i]),N - SA[ip],true}; 
				
      // reset LCP value
      m = std::numeric_limits<int64_t>::max();
      // increment number of runs
      bwtruns++;
		}
	}

  // evaluate last active candidates
  eval(sigma,-1,R,S);

  // remove chached files
  sdsl::remove(cache_file_name(conf::KEY_TEXT, cc));
  sdsl::remove(cache_file_name(conf::KEY_SA, cc));
  sdsl::remove(cache_file_name(conf::KEY_ISA, cc));
  sdsl::remove(cache_file_name(conf::KEY_LCP, cc));

  if(sort) std::sort(S.begin(),S.end());

  if(output_file.length()==0){
    for(auto x : S) cout << x << " ";
    cout << endl;
  }
  else
  {
    ofstream ofs(output_file, ios::binary);

    for (size_t i=0;i<sigma;++i)
    {
    	for (uint64_t s : S[i])
        	ofs.write((char*)&s,5);
    }

    ofs.close();
  }

  if(rho)  cout << "Size of smallest suffixient set: " << S.size() << endl;
  if(runs) cout << "Number of equal-letter BWT(rev(T)) runs: " << bwtruns << endl;
}