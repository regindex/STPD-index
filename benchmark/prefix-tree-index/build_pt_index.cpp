#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <variant>

#include "pt-index.hpp"

void help(){

    std::cout << "build_store_stpd_index [options]" << std::endl <<
    "Options:" << std::endl <<
    "-h          Print usage info." << std::endl <<
    "-i <arg>    Input text file path. (REQUIRED)" << std::endl <<
    "-a <arg>    Prefix array file path (5 bytes per entry). (REQUIRED)" << std::endl <<
    "-s <arg>    LCS-array file path (5 bytes per entry). (REQUIRED)" << std::endl <<
    "-t <arg>    Text oracle variant: (RLZ|bitpacked).         (Def. RLZ)" << std::endl <<
    "-l <arg>    RLZ reference sequence length (if known). (Def. None)" << std::endl <<
    "-o <arg>    Output index file path. (REQUIRED)" << std::endl;
    exit(0);
} 

// index variants
using indexVariant = std::variant<
    stpd::pt_index<RLZ_DNA_sux<>>,
    stpd::pt_index<stpd::bitpacked_text_oracle>
>;

int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        std::cerr << "Wrong number of parameters... See the help messagge:" << std::endl;
        help();
        exit(1);
    }

    std::string inputPath, paPath, lcsPath, outputPath, oracleVariant = "RLZ";
    bool verbose = false;
    size_t refLen = 0;

    indexVariant index;

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
            case 't':
                oracleVariant = std::string(optarg);
            break;
            default:
                help();
            return -1;
        }
    }

    if(inputPath == "" or paPath == "" or outputPath == ""){ help(); }

    // set the index variant
    if(oracleVariant      == "RLZ")      { index.emplace<0>(); }
    else if(oracleVariant == "bitpacked"){ index.emplace<1>(); }
    else                                 { help(); }    

    std::cout << "\n[INFO] Constructing and storing the Prefix tree index" 
              << " for " << inputPath << "\n" << std::endl;

    {   // select the correct index variant
        std::visit([&](auto& idx) {
            // compute the index
            idx.build(inputPath,paPath,lcsPath,refLen);
            // store the index
            idx.store(outputPath); }, index);
    }

    {   // write index parameters
        std::ofstream out(outputPath + ".PT_param");
        out << oracleVariant << "\n" << refLen << "\n" << inputPath << "\n" << outputPath;
        out.close();
    }

    return 0;
}