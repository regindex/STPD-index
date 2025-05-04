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
    "-i <arg>    Input files base path. (REQUIRED)" << std::endl <<
    "-v <arg>    Index variant: (small|large). (Def. small)" << std::endl <<
    "-o <arg>    Output index filename. (REQUIRED)" << std::endl;
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

    std::string inputPath, outputPath, indexVariant;
    bool verbose = false;

    int opt;
    while ((opt = getopt(argc, argv, "hi:o:v:")) != -1)
    {
        switch (opt){
            case 'h':
                help();
            break;
            case 'i':
                inputPath = std::string(optarg);
            break;
            case 'o':
                outputPath = std::string(optarg);
            break;
            case 'v':
                indexVariant = std::string(optarg);
            break;
            default:
                help();
            return -1;
        }
    }

    if(inputPath == "" or outputPath == ""){ help(); }
    
    if(indexVariant == "small")
    {
        std::cout << "### Constructing the ST colex- index for " 
                  << inputPath << std::endl;

        stpd::stpd_index<stpd::stpd_array_binary_search<RLZ_DNA<>>,
                         RLZ_DNA<>,stpd::r_index_phi_inv<>> index;
        index.build_colex_m(inputPath);
        index.store(outputPath);
    }
    else if(indexVariant == "large")
    {
        std::cout << "### Constructing the ST colex+- index for " 
                  << inputPath << std::endl;

        stpd::stpd_index<stpd::stpd_array_binary_search<RLZ_DNA<>>,
                         RLZ_DNA<>,stpd::r_index_phi_inv<>> index;
        index.build_colex_pm(inputPath);
        index.store(outputPath);
    }
    else{ std::cerr << "Not yet implemented..." << std::endl; exit(1); }

    /*
    std::cout << "##############################################" << std::endl;
    std::cout << "##############################################" << std::endl;

    stpd::stpd_index<stpd::stpd_array_binary_search<stpd::bitpacked_text_oracle>,
                     stpd::bitpacked_text_oracle,stpd::stpd_index_phi_inv<>> index__;
    index__.build_colex_m_v2(inputPath);
    index__.store(inputPath+"colex+_v2");
    */

    return 0;
}