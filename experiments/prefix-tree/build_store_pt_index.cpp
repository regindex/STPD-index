#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <unistd.h>

#include "pt-index.hpp"

void help(){

    std::cout << "build_store_stpd_index [options]" << std::endl <<
    "Options:" << std::endl <<
    "-h          Print usage info." << std::endl <<
    "-i <arg>    Input text file path. (REQUIRED)" << std::endl <<
    "-a <arg>    Prefix array file path (5 bytes per entry). (REQUIRED)" << std::endl <<
    "-s <arg>    LCS-array file path (5 bytes per entry). (REQUIRED)" << std::endl <<
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

    std::string inputPath, paPath, lcsPath, outputPath;
    bool verbose = false;
    size_t refLen = 0;

    int opt;
    while ((opt = getopt(argc, argv, "hi:o:l:a:s:")) != -1)
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
            case 's':
                lcsPath = std::string(optarg);
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

    std::cout << "\n[INFO] Constructing and storing the Prefix tree index" 
              << " for " << inputPath << "\n" << std::endl;

    stpd::pt_index<RLZ_DNA_sux<>> index;
    index.build(inputPath,paPath,lcsPath,refLen);
    // store the index
    index.store(outputPath);

    return 0;
}