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
    "-t <arg>    Maximum number of occurrences to report per pattern. (Def. all)" << std::endl <<
    "-e <arg>    Locate all exponential search base. (Def. 2)" << std::endl <<
    "-c          Check occurrences correctness. (Def. False)" << std::endl <<
    "-b          Run locate queries performance benchmark (Def. False)" << std::endl;
    exit(0);
} 

// index variants
using indexVariant = std::variant<
    stpd::colex_minus_index<stpd::tabulated_binary_search_DNA<RLZ_DNA_sux<>>,
                            RLZ_DNA_sux<>,
                            stpd::r_index_phi_inv>,
    stpd::colex_minus_index<stpd::tabulated_binary_search_DNA<stpd::bitpacked_text_oracle>,
                            stpd::bitpacked_text_oracle,
                            stpd::r_index_phi_inv>,
    stpd::colex_minus_index<stpd::tabulated_binary_search_DNA<RLZ_DNA_sux<>>,
                            RLZ_DNA_sux<>,
                            stpd::move_r_phi_inv>,
    stpd::colex_minus_index<stpd::tabulated_binary_search_DNA<stpd::bitpacked_text_oracle>,
                            stpd::bitpacked_text_oracle,
                            stpd::move_r_phi_inv>
>;

int main(int argc, char* argv[])
{
    std::setlocale(LC_ALL, "C");

    if(argc < 2)
    {
        std::cerr << "Wrong number of parameters... See the help messagge:" << std::endl;
        help();
        exit(1);
    }

    std::string input_path, pattern_file;
    usafe_t max_occ_thr = std::numeric_limits<usafe_t>::max();
    double exp_search_base = 2;
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
            case 't':
                max_occ_thr = std::stoull(optarg);
            break;
            case 'e':
                exp_search_base = std::stod(optarg);
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
        std::string phiVariant, oracleVariant;
        std::ifstream in(input_path + ".param");
        in >> phiVariant;
        in >> oracleVariant;
        in.close();

        // set the index variant
        if(oracleVariant      == "RLZ"       and phiVariant == "r-index"){ index.emplace<0>(); }
        else if(oracleVariant == "bitpacked" and phiVariant == "r-index"){ index.emplace<1>(); }
        else if(oracleVariant == "RLZ"       and phiVariant == "move")   { index.emplace<2>(); }
        else if(oracleVariant == "bitpacked" and phiVariant == "move")   { index.emplace<3>(); }
        else                                                             { help(); }
    }

    // Load STPD-index from file
    std::visit([&](auto& idx) {
        // load the index
        if( max_occ_thr == 1 )
            idx.load(input_path, true);
        else
            idx.load(input_path);
    }, index);

    // Run locate queries benchmark
    if( run_benchmark )
    {
        if( max_occ_thr == 1 ) {
            std::cout << "### Running locate primary occurrence queries benchmark for "
                      << pattern_file << " using the index in "
                      << input_path << std::endl;   

            std::visit([&](auto& idx) {
                    idx.locate_primary_benchmark(pattern_file); 
            }, index); 
        }
        else {
            std::cout << "### Running locate all occurrence queries benchmark for "
                      << pattern_file << " using the index in "
                      << input_path << std::endl;

            std::visit([&](auto& idx) {
                    idx.locate_secondary_benchmark(pattern_file, max_occ_thr, exp_search_base); 
            }, index);
        }

        return 0;
    }

    std::cout << "### Running locate queries for " << pattern_file 
              << " using the index in " << input_path << std::endl;

    std::visit([&](auto& idx) {
        if( max_occ_thr == 1 )
            idx.locate_one_fasta(pattern_file, check_correctness); 
        else
            idx.locate_all_fasta(pattern_file, max_occ_thr, exp_search_base, check_correctness); 
    }, index);

    return 0;
}