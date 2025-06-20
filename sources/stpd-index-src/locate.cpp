#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <unistd.h>

#include "stpd-index.hpp"

void help(){

    std::cout << "locate [options]" << std::endl <<
    "Options:" << std::endl <<
    "-h          Print usage info." << std::endl <<
    "-i <arg>    Input index filepath. (REQUIRED)" << std::endl <<
    "-p <arg>    Patterns FASTA file.  (REQUIRED)" << std::endl <<
    "-t <arg>    Maximum number of occurrences to report per pattern. (Def. none)" << std::endl;
    //"-O <arg>    Enable DNA index optimizations: (v1|v2|v3). (Def. False)" << std::endl;
    exit(0);
} 

int main(int argc, char* argv[])
{
    if(argc < 3)
    {
        std::cerr << "Wrong number of parameters... See the help messagge:" << std::endl;
        help();
        exit(1);
    }

    std::string inputPath, patternFile; //optVariant;
    uint64_t maxOcc = (1ULL << 63) | ((1ULL << 63) - 1);
    bool verbose = false;

    int opt;
    while ((opt = getopt(argc, argv, "hi:p:O:t:")) != -1)
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
            case 't':
                maxOcc = std::stoull(optarg);
            break;
            //case 'O':
            //    optVariant = std::string(optarg);
            //break;
            default:
                help();
            return -1;
        }
    }

    {
        std::cout << "### Running locate all occurrence queries for "
                  << patternFile << " using the index in "
                  << inputPath << std::endl;

        stpd::stpd_index<stpd::stpd_array_binary_search_opt<>,
                         RLZ_DNA_sux<>,stpd::r_index_phi_inv_intlv> index;
        index.load(inputPath);
        // run locate all occurrence queries
        index.locate_fasta(patternFile,maxOcc);
    }

    return 0;
}