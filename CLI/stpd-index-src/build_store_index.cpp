#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <filesystem>
#include <variant>

#include <colex-index.hpp>

void help(){

    std::cout << "stpd-build [options]" << std::endl <<
    "Options:" << std::endl <<
    "-i <arg>    Input text base path.                   (REQUIRED)" << std::endl <<
    "-o <arg>    Output index file path.                 (REQUIRED)" << std::endl <<
    "-p <arg>    Phi-function machinery: (r-index|move). (Def. r-index)" << std::endl <<
    "-t <arg>    Text oracle variant: (RLZ|bitpacked).   (Def. RLZ)" << std::endl <<
    "-b <arg>    Tabulation sequences length.            (Def. 15)" << std::endl <<
    "-l <arg>    RLZ reference prefix length.            (Def. None)" << std::endl <<
    "-k          Keep temporary files.                   (Def. False)" << std::endl <<
    "-h          Print usage info." << std::endl;
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

    std::string inputPath, outputPath, 
                phiVariant = "r-index", oracleVariant = "RLZ";
    bool verbose = false, keep = false;
    usafe_t refLen = 0, tabLen = 15;

    indexVariant index;

    int opt;
    while ((opt = getopt(argc, argv, "hi:o:t:p:l:b:k")) != -1)
    {
        switch (opt){
            case 'h':
                help();
            break;
            case 'k':
                keep = true;
            break;
            case 'i':
                inputPath = std::string(optarg);
            break;
            case 'o':
                outputPath = std::string(optarg);
            break;
            case 'p':
                phiVariant = std::string(optarg);
            break;
            case 't':
                oracleVariant = std::string(optarg);
            break;
            case 'l':
                refLen = std::atoll(optarg);
            break;
            case 'b':
                tabLen = std::atoll(optarg);
            break;
            default:
                help();
            return -1;
        }
    }

    // the input and output paths cannot be emtpy
    if(inputPath == "" or outputPath == ""){ help(); }

    // set the index variant
    if(oracleVariant      == "RLZ"   and phiVariant == "r-index")    { index.emplace<0>(); }
    else if(oracleVariant == "bitpacked" and phiVariant == "r-index"){ index.emplace<1>(); }
    else if(oracleVariant == "RLZ"   and phiVariant == "move")       { index.emplace<2>(); }
    else if(oracleVariant == "bitpacked" and phiVariant == "move")   { index.emplace<3>(); }
    else                                                             { help(); }

    if(not std::filesystem::exists(inputPath))
    {
        std::cerr << "Input file: " << inputPath << " not found! exiting..." << std::endl;
        exit(1);
    }

    std::string path_decomposition = inputPath + ".colex_m";
    if( not std::filesystem::exists(path_decomposition) )
    {
        std::cout << "\n[INFO] Constructing and storing the Suffix Tree path decomposition index (STDP-index)" 
                  << " for " << inputPath << " using the " << phiVariant << " phi-function machinery\n" << std::endl;

        { // compute the path decomposition
            std::cout << "[STEP 0] Computing the ST path decomposition..." << "\n" << std::endl;

            std::filesystem::path exe_path = std::filesystem::canonical(argv[0]);
            std::filesystem::path exe_dir = exe_path.parent_path();

            std::string command = std::string(exe_dir) +
                                  "/stpd-sampling -i " + inputPath + " -o " + inputPath + ".colex_m -c -P"; 
            int result = std::system(command.c_str());
            if (result != 0) {
                std::cerr << "Error while computing the path decomposition..." << std::endl;
                exit(1);
            }
        }
    }
    else
        std::cout << "\n[INFO] Suffix Tree path decomposition found! Proceeding to index construction\n" << std::endl;

    { // select the correct index variant
        std::visit([&](auto& idx) {
            // compute the index
            idx.build(inputPath,inputPath+".colex_m",inputPath+".rbwt",
                      inputPath+".pa",inputPath+".lcs",refLen,tabLen);
            // store the index
            idx.store(outputPath); }, index);
    }

    { // write index parameters
        std::ofstream out(outputPath + ".param");
        out << phiVariant << "\n" << oracleVariant << "\n" << refLen << "\n" << inputPath << "\n" << outputPath;
        out.close();
    }

    if( not keep )
    { // delete temporary files
        std::remove(std::string(inputPath+".colex_m").c_str());
        std::remove(std::string(inputPath+".rbwt").c_str());
        std::remove(std::string(inputPath+".pa").c_str());
        std::remove(std::string(inputPath+".lcs").c_str());
        std::remove(std::string(inputPath+".rlz").c_str());
    }

    return 0;
}