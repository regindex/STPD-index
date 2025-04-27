// Copyright (c) 2025, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#include <iostream>
#include <unistd.h>

#include "lcp.hpp"
#include "dna_bwt.hpp"

using namespace std;

string input_bwt, input_sa; 
string output_file;

bool out_da = false;
uint8_t lcp_size = 1;

bool DNAonly = false;

char TERM = '#';

void help(){

    cout << "sampling [options]" << endl <<
    "Input: BWT and SA of a text. Output: ST path decomposition of the text." << endl <<
    "Options:" << endl <<
    "-h          Print this help" << endl <<
    "-i <arg>    Input BWT (REQUIRED)" << endl <<
    "-s <arg>    Input SA (REQUIRED)" << endl <<
    "-o <arg>    Output file name (REQUIRED)" << endl <<
    //"-l <arg>    Number of Bytes used to represent LCP values. <arg>=1,2,4,8 Bytes. Default: 1." << endl <<
    //"-n          Alphabet is {A,C,G,N,T," << TERM << "}. Default: alphabet is {A,C,G,T," << TERM << "}."<< endl <<
    "-t          ASCII code of the terminator. Default:" << int('#') << " (#). Cannot be the code for A,C,G,T,N." << endl;
    exit(0);
}

int main(int argc, char** argv){

    if(argc < 3) help();

    int opt;
    while ((opt = getopt(argc, argv, "hi:o:s:t:")) != -1){
        switch (opt){
            case 'h':
                help();
            break;
            case 'i':
                input_bwt = string(optarg);
            break;
            case 's':
                input_sa  = string(optarg);
            break;
            case 'o':
                output_file = string(optarg);
            break;
            case 't':
                TERM = atoi(optarg);
            break;
            default:
                help();
            return -1;
        }
    }

    if(TERM == 'A' or TERM == 'C' or TERM == 'G' or TERM == 'T')
    {
        cout << "Error: invalid terminator '" << TERM << "'" << endl;
        help();
    }
    if(input_bwt.size()==0) help();
    if(output_file.size()==0) help();

    cout << "Input BWT file: " << input_bwt << endl;
    cout << "Input SA file: " << input_sa << endl;
    cout << "Output sampling file: " << output_file << endl;

    DNAonly = hasDNAonly(input_bwt,TERM);
    std::cout << int(DNAonly) << std::endl;

    if(DNAonly){

        cout << "Alphabet: A,C,G,T,'" << TERM << "'" << endl;

        cout << "Loading and indexing BWT ... " << endl;

        dna_bwt_t BWT = dna_bwt_t(input_bwt, TERM);

        cout << "Done. Size of BWT: " << BWT.size() << endl;

        lcp<dna_bwt_t, uint64_t> M1(&BWT);
        cout << "Storing output to file ... " << endl;
        M1.save_to_file(output_file);

    }

    cout << "Done. " << endl;

}