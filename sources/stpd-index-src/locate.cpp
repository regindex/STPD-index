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
    "-p <arg>    Patterns FASTA file.  (REQUIRED)" << std::endl <<
    "-O <arg>    Enable DNA index optimizations: (v1|v2). (Def. False)" << std::endl;
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

    std::string inputPath, patternFile, optVariant;
    bool verbose = false;

    int opt;
    while ((opt = getopt(argc, argv, "hi:p:O:")) != -1)
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
            case 'O':
                optVariant = std::string(optarg);
            break;
            default:
                help();
            return -1;
        }
    }
    
    if(optVariant == "v1")
    {
        std::cout << "### Querying DNA optimized ST colex index for " 
                  << inputPath << std::endl;

        stpd::stpd_index<stpd::stpd_array_binary_search_opt<RLZ_DNA<>>,
                         RLZ_DNA<>,stpd::r_index_phi_inv<>> index;
        index.load(inputPath);

        index.locate_fasta(patternFile);
    }
    else if(optVariant == "v2")
    {
        std::cout << "### Querying DNA optimized v2 ST colex index for " 
                  << inputPath << std::endl;

        stpd::stpd_index<stpd::stpd_array_binary_search_opt_v2<RLZ_DNA<>>,
                         RLZ_DNA<>,stpd::r_index_phi_inv_sux> index;
        index.load(inputPath);

        index.locate_fasta(patternFile);
    }
    else
    {
        std::cout << "### Querying ST colex index for " 
                  << inputPath << std::endl;

        stpd::stpd_index<stpd::stpd_array_binary_search<RLZ_DNA<>>,
                         RLZ_DNA<>,stpd::r_index_phi_inv<>> index;
        index.load(inputPath);

        index.locate_fasta(patternFile);
    }

    return 0;
}