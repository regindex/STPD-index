#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <filesystem>
#include <variant>

#include "pa-index.hpp"

void help(){

    std::cout << "locate [options]" << std::endl <<
    "Options:" << std::endl <<
    "-h          Print usage info." << std::endl <<
    "-i <arg>    Input index filepath. (REQUIRED)" << std::endl <<
    "-p <arg>    Patterns FASTA file.  (REQUIRED)" << std::endl <<
    "-t <arg>    Maximum number of occurrences to report per pattern. (Def. all)" << std::endl <<
    "-b          Run queries in benchmark mode. (Def. false)" << std::endl;
    exit(0);
} 

// index variants
using indexVariant = std::variant<
    stpd::pa_index<RLZ_DNA_sux<>>,
    stpd::pa_index<stpd::bitpacked_text_oracle>
>;

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
    bool benchmark = false;

    indexVariant index;

    int opt;
    while ((opt = getopt(argc, argv, "hi:p:t:b")) != -1)
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
            case 'b':
                benchmark = true;
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

    { // read index parameters
        std::string oracleVariant;
        std::ifstream in(inputPath + ".PA_param");
        in >> oracleVariant;
        in.close();

        // set the index variant
        if(oracleVariant      == "RLZ")  { index.emplace<0>(); }
        else if(oracleVariant == "plain"){ index.emplace<1>(); }
        else                             { help(); }   
    }

    // Load PA-index from file
    std::visit([&](auto& idx) {
        // load the index
        idx.load(inputPath); }, index);

    std::cout << "### Running locate all occurrence queries for "
              << patternFile << " using the index in "
              << inputPath << std::endl;

    // run locate queries
    std::visit([&](auto& idx) {
        if( benchmark )
            idx.locate_fasta_benchmark(patternFile,maxOcc); 
        else
            idx.locate_fasta(patternFile,maxOcc); 
    }, index);

    return 0;
}