#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <unistd.h>

#include "stpd-index.hpp"

void help(){

    std::cout << "build_DNA_index [options]" << std::endl <<
    "Options:" << std::endl <<
    "-h          Print usage info." << std::endl <<
    "-i <arg>    Input index filepath. (REQUIRED)" << std::endl <<
    "-p <arg>    Patterns FASTA file.  (REQUIRED)" << std::endl;
    exit(0);
} 

int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        std::cerr << "Wrong number of parameters... See the help messagge:" << std::endl;
        help();
        exit(1);
    }

    std::string inputPath, patternFile;
    bool verbose = false;

    int opt;
    while ((opt = getopt(argc, argv, "hi:p:")) != -1)
    {
        switch (opt){
            case 'h':
                help();
            break;
            case 'i':
                inputPath = std::string(optarg);
            break;
            case 'p':
                patternFile = std::string(optarg);
            break;
            default:
                help();
            return -1;
        }
    }
    
    std::cout << "### Querying ST colex index for " 
              << inputPath << std::endl;

    stpd::stpd_index<stpd::stpd_array_binary_search<RLZ_DNA<>>,
                     RLZ_DNA<>,stpd::r_index_phi_inv<>> index;
    index.load(inputPath);

    index.locate_fasta(patternFile);

    return 0;
}