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
    "-v <arg>    Index variant: (colex-|colex+-). (Def. colex-)" << std::endl <<
    "-O          Enable index optimizations. (Def. False)" << std::endl <<
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
    bool verbose = false, optimizations = false;

    int opt;
    while ((opt = getopt(argc, argv, "hi:o:v:O")) != -1)
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
            case 'O':
                optimizations = true;
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
    
    if(indexVariant == "colex-" and not optimizations)
    {
        std::cout << "### Constructing the ST colex- index for " 
                  << inputPath << std::endl;

        stpd::stpd_index<stpd::stpd_array_binary_search<RLZ_DNA<>>,
                         RLZ_DNA<>,stpd::r_index_phi_inv<>> index;
        index.build_colex_m(inputPath,inputPath+".colex_m",inputPath+".rbwt",
                            inputPath+".pa");
        index.store(outputPath);
    }
    else if(indexVariant == "colex+-" and not optimizations)
    {
        std::cout << "### Constructing the ST colex+- index for " 
                  << inputPath << std::endl;

        stpd::stpd_index<stpd::stpd_array_binary_search<RLZ_DNA<>>,
                         RLZ_DNA<>,stpd::r_index_phi_inv<>> index;
        index.build_colex_pm(inputPath,inputPath+".colex_pm",inputPath+".rbwt",
                             inputPath+".pa");
        index.store(outputPath);
    }
    else if(indexVariant == "colex-" and optimizations)
    {
        std::cout << "### Constructing the opt ST colex- index for " 
                  << inputPath << std::endl;

        stpd::stpd_index<stpd::stpd_array_binary_search_opt<RLZ_DNA<>>,
                         RLZ_DNA<>,stpd::r_index_phi_inv<>> index;
        index.build_colex_m(inputPath,inputPath+".colex_m",inputPath+".rbwt",
                            inputPath+".pa",inputPath+".lcs");
        index.store(outputPath);
    }
    else{ std::cerr << "Not yet implemented..." << std::endl; exit(1); }

    return 0;
}