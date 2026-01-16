#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <filesystem>
#include <variant>

#include <colex-index.hpp>

void help(){

    std::cout << "stpd-locate [options]" << std::endl <<
    "Options:" << std::endl <<
    "-h          Print usage info." << std::endl <<
    "-i <arg>    Input index filepath. (REQUIRED)" << std::endl <<
    "-p <arg>    Patterns FASTA file.  (REQUIRED)" << std::endl <<
    "-t <arg>    Maximum number of occurrences to report per pattern. (Def. none)" << std::endl <<
    "-c          Check occurrences correctness. (Def. False)" << std::endl <<
    "-b          Run queries in benchmark mode. (Def. False)" << std::endl;
    exit(0);
} 

// index variants
using indexVariant = std::variant<
    stpd::colex_minus_index<stpd::tabulated_binary_search_DNA<RLZ_DNA_sux<>>,
                            RLZ_DNA_sux<>,
                            stpd::r_index_phi_inv_intlv>,
    stpd::colex_minus_index<stpd::tabulated_binary_search_DNA<stpd::bitpacked_text_oracle>,
                            stpd::bitpacked_text_oracle,
                            stpd::r_index_phi_inv_intlv>,
    stpd::colex_minus_index<stpd::tabulated_binary_search_DNA<RLZ_DNA_sux<>>,
                            RLZ_DNA_sux<>,
                            stpd::move_r_phi_inv>,
    stpd::colex_minus_index<stpd::tabulated_binary_search_DNA<stpd::bitpacked_text_oracle>,
                            stpd::bitpacked_text_oracle,
                            stpd::move_r_phi_inv>
>;

int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        std::cerr << "Wrong number of parameters... See the help messagge:" << std::endl;
        help();
        exit(1);
    }

    std::string inputPath, patternFile;
    uint64_t maxOcc = (1ULL << 63) | ((1ULL << 63) - 1);
    bool check_correctness = false, run_benchmark = false;

    indexVariant index;

    int opt;
    while ((opt = getopt(argc, argv, "hi:p:t:cb")) != -1)
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

    if(inputPath == "" or patternFile == ""){ help(); }

    if(not std::filesystem::exists(inputPath))
    {
        std::cerr << "Index file: " << inputPath << " not found! exiting..." << std::endl;
        exit(1);
    }

    { // read index parameters
        std::string phiVariant, oracleVariant;
        std::ifstream in(inputPath + ".param");
        in >> phiVariant;
        in >> oracleVariant;
        in.close();

        // set the index variant
        if(oracleVariant      == "RLZ"   and phiVariant == "r-index"){ index.emplace<0>(); }
        else if(oracleVariant == "plain" and phiVariant == "r-index"){ index.emplace<1>(); }
        else if(oracleVariant == "RLZ"   and phiVariant == "move")   { index.emplace<2>(); }
        else if(oracleVariant == "plain" and phiVariant == "move")   { index.emplace<3>(); }
        else                                                         { help(); }
    }

    // Load STPD-index from file
    std::visit([&](auto& idx) {
        // load the index
        idx.load(inputPath); }, index);

    if(maxOcc == 1){
        std::cout << "### Running locate one occurrence queries for "
                  << patternFile << " using the index in "
                  << inputPath << std::endl;

        // run locate one occurrence queries
        std::visit([&](auto& idx) {
            if( run_benchmark )
                idx.locate_one_benchmark(patternFile, check_correctness); 
            else
                idx.locate_one_fasta(patternFile, check_correctness); 
        }, index);
    }
    else{
        std::cout << "### Running locate all occurrence queries benchmark for "
                  << patternFile << " using the index in "
                  << inputPath << std::endl;

        // run locate all occurrence queries
        std::visit([&](auto& idx) {
            if( run_benchmark )
                idx.locate_all_benchmark(patternFile, maxOcc, check_correctness);
            else
                idx.locate_all_fasta(patternFile, maxOcc, check_correctness);
        }, index);
    }

    return 0;
}