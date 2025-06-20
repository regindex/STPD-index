#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <unistd.h>

#include "stpd-index.hpp"

void help(){

    std::cout << "build_store_stpd_index [options]" << std::endl <<
    "Options:" << std::endl <<
    "-h          Print usage info." << std::endl <<
    "-i <arg>    Input text file path. (REQUIRED)" << std::endl <<
    //"-v <arg>    Index variant: (colex-|colex+-). (REQUIRED)" << std::endl <<
    //"-O <arg>    Enable DNA index optimizations: (v1|v2|v3). (Def. False)" << std::endl <<
    "-l <arg>    RLZ reference sequence length (if known). (Def. None)" << std::endl <<
    "-o <arg>    Output index file path. (REQUIRED)" << std::endl;
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

    std::string inputPath, outputPath; // indexVariant, optVariant;
    bool verbose = false;
    size_t refLen = 0;

    int opt;
    while ((opt = getopt(argc, argv, "hi:o:v:O:l:")) != -1)
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
            //case 'O':
            //    optVariant = std::string(optarg);
            //break;
            //case 'v':
            //    indexVariant = std::string(optarg);
            //break;
            case 'l':
                refLen = std::atoll(optarg);
            break;
            default:
                help();
            return -1;
        }
    }

    if(inputPath == "" or outputPath == ""){ help(); }

    std::cout << "### Constructing and storing the Suffix Tree path decomposition (STDP)" 
              << " index for " << inputPath << std::endl;

    { // compute the path decomposition
        std::string command = "./sources/path-decomp-src/stpd_small -i " + inputPath + " -o " + inputPath + ".colex_m -c -P"; 
        int result = std::system(command.c_str());
        if (result != 0) {
            std::cerr << "Error while computing the path decomposition..." << std::endl;
            exit(1);
        }
    }

    { // compute the index
        stpd::stpd_index<stpd::stpd_array_binary_search_opt<>,
                         RLZ_DNA_sux<>,stpd::r_index_phi_inv_intlv> index;
        index.build_colex_m(inputPath,inputPath+".colex_m",inputPath+".rbwt",
                            inputPath+".pa",inputPath+".lcs",refLen);
        // store the index
        index.store(outputPath);
    }

    { // delete temporary files
        std::remove(std::string(inputPath+".colex_m").c_str());
        std::remove(std::string(inputPath+".rbwt").c_str());
        std::remove(std::string(inputPath+".pa").c_str());
        std::remove(std::string(inputPath+".lcs").c_str());
        std::remove(std::string(inputPath+".rlz").c_str());
    }

    return 0;
}