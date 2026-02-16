#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <filesystem>
#include <variant>

#include "sA-index.hpp"

void help(){

    std::cout << "stpd-locate [options]" << std::endl <<
    "Options:" << std::endl <<
    "-h          Print usage info." << std::endl <<
    "-i <arg>    Input index filepath. (REQUIRED)" << std::endl <<
    "-p <arg>    Patterns FASTA file.  (REQUIRED)" << std::endl <<
    "-c          Check occurrences correctness. (Def. False)" << std::endl <<
    "-b          Run locate queries performance benchmark (Def. False)" << std::endl;
    exit(0);
} 

// index variants
using indexVariant = std::variant<
    suffixient::suffixient_array_index<stpd::tabulated_binary_search_DNA<RLZ_DNA_sux<>>,
                                       RLZ_DNA_sux<>>,
    suffixient::suffixient_array_index<stpd::tabulated_binary_search_DNA<stpd::bitpacked_text_oracle>,
                                       stpd::bitpacked_text_oracle>
>;

int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        std::cerr << "Wrong number of parameters... See the help messagge:" << std::endl;
        help();
        exit(1);
    }

    std::string input_path, pattern_file;
    bool check_correctness = false, run_benchmark = false;

    indexVariant index;

    int opt;
    while ((opt = getopt(argc, argv, "hi:p:t:e:cb")) != -1)
    {
        switch (opt){
            case 'h':
                help();
            break;
            case 'i':
                input_path = std::string(optarg);
            break;
            case 'p':
                pattern_file = std::string(optarg);
            break;
            case 'c':
                check_correctness = true;
            break;
            case 'b':
                run_benchmark = true;
            break;
            default:
                help();
            return -1;
        }
    }

    if(input_path == "" or pattern_file == ""){ help(); }

    if(not std::filesystem::exists(input_path))
    {
        std::cerr << "Index file: " << input_path << " not found! exiting..." << std::endl;
        exit(1);
    }

    { // read index parameters
        std::string oracleVariant;
        std::ifstream in(input_path + ".param");
        in >> oracleVariant;
        in.close();
        std::cout << oracleVariant << std::endl;

        // set the index variant
        if(oracleVariant      == "RLZ")       { index.emplace<0>(); }
        else if(oracleVariant == "bitpacked") { index.emplace<1>(); }
        else                                  { help(); }
    }

    // Load suffixient-array index from file
    std::visit([&](auto& idx) {
        // load the index
        idx.load(input_path);
    }, index);

    // Run locate queries benchmark
    if( run_benchmark )
    {
        std::cout << "### Running locate primary occurrence queries benchmark for "
                  << pattern_file << " using the index in "
                  << input_path << std::endl;   

        std::visit([&](auto& idx) {
                idx.locate_primary_benchmark(pattern_file); 
        }, index); 

        return 0;
    }
    /*
    std::cout << "### Running locate queries for " << pattern_file 
              << " using the index in " << input_path << std::endl;

    std::visit([&](auto& idx) {
        if( max_occ_thr == 1 )
            idx.locate_one_fasta(pattern_file, check_correctness); 
        else
            idx.locate_all_fasta(pattern_file, max_occ_thr, exp_search_base, check_correctness); 
    }, index);
    */

    return 0;
}