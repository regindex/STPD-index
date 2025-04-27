// Copyright (c) 2025, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#include <iostream>
#include <sdsl/construct.hpp>
#include <set>
#include <limits>
#include <algorithm>

using namespace std;
using namespace sdsl;

#define STORE_SIZE 5

int_vector<8> T;
int_vector<40> A;
int_vector_buffer<> SA;

cache_config cc;
uint64_t N = 0; //including 0x0 terminator

inline uint8_t BWT(uint64_t i){ return SA[i] == 0 ? 0 : T[SA[i] - 1]; }

void help()
{
	cout << "Simple software to compute the PA and BWT of a text [options]" << endl <<
	"Input: non-empty ASCII file without character 0x0, from standard input." << endl <<
	"Output: The PA of the input and the BWT of the reversed input." << endl <<
	"Warning: if 0x0 appears, the standard input is read only until the first occurrence of 0x0 (excluded)." << endl <<
	"Options:" << endl <<
	"-h          Print usage info." << endl <<
	//"-b          Compute the BWT and SA of the input text" << endl <<
	//"-r          Compute the BWT and SA (PA) of the reversed input text" << endl <<
	//"-s          Compute the Suffix Array of the input text" << endl <<
	//"-p          Compute the Prefix Array of the input text" << endl <<
	"-o <arg>    Output files basepath." << endl;
	exit(0);
} 

int main(int argc, char** argv)
{
	if(argc < 2) help();

	string output_basepath;
	bool bwt = false, bwt_r = false, sa = false, pa = false;

	int opt;
	while ((opt = getopt(argc, argv, "ho:")) != -1)
	{
		switch (opt)
		{
			case 'h':
				help();
			break;
			case 'o':
				output_basepath = string(optarg);
			break;
			default:
				help();
			return -1;
		}
	}

	{
		string in;
		getline(cin,in,char(0));
		N = in.size() + 1;

		if(N<2){
			cerr << "Error: empty text" <<  endl;
			help();
		}

		T = int_vector<8>(N - 1);
		A = int_vector<40>(128,0);

		{
			for(uint64_t i = 0; i < N - 1; ++i)
			{
				uint8_t c = in[N - i - 2];
				T[i] = c;
				A[c]++;
			}
		}

		append_zero_symbol(T);
		store_to_cache(T, conf::KEY_TEXT, cc);
		construct_sa<8>(cc);
		SA = int_vector_buffer<>(cache_file_name(conf::KEY_SA, cc));
	}

	{
	  	ofstream ofs = ofstream(output_basepath + ".pa", ios::binary);
	  	for (size_t i=0;i<N;++i)
	  		{ size_t x = N - SA[i] - 1; ofs.write(reinterpret_cast<const char*>(&x), STORE_SIZE); }
	  	ofs.close();
	  	ofs = ofstream(output_basepath + ".rbwt", ios::binary);
	  	for (size_t i=0;i<N;++i)
	  		{ char x = T[SA[i]-1]; ofs.write(&x, 1); }
	  	ofs.close();
	  	ofs = ofstream(output_basepath + ".mult", ios::binary);
	  	for (size_t i=0;i<128;++i)
	  		{ size_t x = A[i]; ofs.write(reinterpret_cast<const char*>(&x), STORE_SIZE); }
	  	ofs.close();
  	}

	// remove chached files
	sdsl::remove(cache_file_name(conf::KEY_TEXT,cc));
	sdsl::remove(cache_file_name(conf::KEY_SA,  cc));
	
	return 0;
}