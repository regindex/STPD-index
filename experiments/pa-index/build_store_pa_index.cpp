#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <unistd.h>

#include "pa-index.hpp"

void help(){

    std::cout << "build_store_stpd_index [options]" << std::endl <<
    "Options:" << std::endl <<
    "-h          Print usage info." << std::endl <<
    "-i <arg>    Input text file path. (REQUIRED)" << std::endl <<
    "-a <arg>    Prefix array file path (5 bytes per entry). (REQUIRED)" << std::endl <<
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

    std::string inputPath, paPath, outputPath;
    bool verbose = false;
    size_t refLen = 0;

    int opt;
    while ((opt = getopt(argc, argv, "hi:o:l:a:")) != -1)
    {
        switch (opt){
            case 'h':
                help();
            break;
            case 'i':
                inputPath = std::string(optarg);
            break; 
            case 'a':
                paPath = std::string(optarg);
            break;
            case 'o':
                outputPath = std::string(optarg);
            break;
            case 'l':
                refLen = std::atoll(optarg);
            break;
            default:
                help();
            return -1;
        }
    }

    if(inputPath == "" or paPath == "" or outputPath == ""){ help(); }

    std::cout << "\n[INFO] Constructing and storing the Suffix Array index" 
              << " for " << inputPath << "\n" << std::endl;

    stpd::pa_index<RLZ_DNA_sux<>> index;
    index.build(inputPath,paPath,refLen);
    // store the index
    index.store(outputPath);

    /*

    { // compute the path decomposition
        std::cout << "[STEP 0] Computing the ST path decomposition..." << "\n" << std::endl;
        std::string command = "./build/sources/path-decomp-src/stpd_small -i " + inputPath + " -o " + inputPath + ".colex_m -c -P"; 
        int result = std::system(command.c_str());
        if (result != 0) {
            std::cerr << "Error while computing the path decomposition..." << std::endl;
            exit(1);
        }
    }

    if(phiDS == "r-index")
    { // compute the index
        stpd::stpd_index<stpd::stpd_array_binary_search_opt<>,
                         RLZ_DNA_sux<>,stpd::r_index_phi_inv_intlv> index;
        index.build_colex_m(inputPath,inputPath+".colex_m",inputPath+".rbwt",
                            inputPath+".pa",inputPath+".lcs",refLen);
        // store the index
        index.store(outputPath);
    }
    else if(phiDS == "move")
    {
        std::cout << "Build move!" << std::endl;
        stpd::stpd_index<stpd::stpd_array_binary_search_opt<>,
                         RLZ_DNA_sux<>,stpd::move_r_phi_inv> index;
        index.build_colex_m(inputPath,inputPath+".colex_m",inputPath+".rbwt",
                            inputPath+".pa",inputPath+".lcs",refLen);
        // store the index
        index.store(outputPath);
    }
    else{ std::cerr << "Please, select a valid phi-function data structure!" << std::endl; exit(1); }

    { // delete temporary files
        std::remove(std::string(inputPath+".colex_m").c_str());
        std::remove(std::string(inputPath+".rbwt").c_str());
        std::remove(std::string(inputPath+".pa").c_str());
        std::remove(std::string(inputPath+".lcs").c_str());
        std::remove(std::string(inputPath+".rlz").c_str());
    }
    */

    return 0;
}