#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <filesystem>

#include "stpd-index.hpp"

void help(){

    std::cout << "locate [options]" << std::endl <<
    "Options:" << std::endl <<
    "-h          Print usage info." << std::endl <<
    "-i <arg>    Input index filepath. (REQUIRED)" << std::endl <<
    "-p <arg>    Patterns FASTA file.  (REQUIRED)" << std::endl <<
    "-t <arg>    Maximum number of occurrences to report per pattern. (Def. none)" << std::endl <<
    "-B          Run locate all queries benchmark. (Def. False, doesn't generate output files)" << std::endl;
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

    std::string inputPath, patternFile;
    uint64_t maxOcc = (1ULL << 63) | ((1ULL << 63) - 1);
    bool verbose = false, bench = false;

    int opt;
    while ((opt = getopt(argc, argv, "hi:p:Bt:")) != -1)
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
            case 'B':
                bench = true;
            break;
            default:
                help();
            return -1;
        }
    }

    if(inputPath == "" or patternFile == ""){ help(); }

    if(not std::filesystem::exists(inputPath))
    {
        std::cerr << "Index file: " << inputPath << " not found! exiting..." << std::endl;
        exit(1);
    }

    // Load STPD-index from file
    stpd::stpd_index<stpd::stpd_array_binary_search_opt<>,
                     RLZ_DNA_sux<>,stpd::r_index_phi_inv_intlv> index;
    index.load(inputPath);

    if(not bench){
        std::cout << "### Running locate all occurrence queries for "
                  << patternFile << " using the index in "
                  << inputPath << std::endl;

        // run locate all occurrence queries
        index.locate_fasta(patternFile,maxOcc);
    }
    else{
        std::cout << "### Running locate all occurrence queries benchmark for "
                  << patternFile << " using the index in "
                  << inputPath << std::endl;

        // run locate all occurrence queries benchmark
        index.locate_fasta_test_running_time(patternFile);
    }

    return 0;
}