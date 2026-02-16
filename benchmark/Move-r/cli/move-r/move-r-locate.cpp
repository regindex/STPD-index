/**
 * part of LukasNalbach/Move-r
 *
 * MIT License
 *
 * Copyright (c) Lukas Nalbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/* Modified by Davide Cenzato to implement the STPD-index benchmarks */

#include <filesystem>
#include <iostream>
#include <move_r/move_r.hpp>

int ptr = 1;
bool output_occurrences = false;
bool check_correctness = false;
uint64_t max_occ = (1ULL << 63) | ((1ULL << 63) - 1);
std::string input;
std::ofstream mf;
std::string path_index_file;
std::string path_patterns_file;
std::string path_output_file;
std::ifstream index_file;
std::ifstream patterns_file;
std::ifstream input_file;
std::ofstream output_file;
std::string name_text_file;
std::string path_input_file;

void help(std::string msg)
{
    if (msg != "") std::cout << msg << std::endl;
    std::cout << "move-r-locate: locate all occurrences of the input patterns." << std::endl << std::endl;
    std::cout << "usage: move-r-locate [options] <index_file> <patterns>" << std::endl;
    std::cout << "   -m <m_file> <text_name>    m_file is the file to write measurement data to," << std::endl;
    std::cout << "                              text_name should be the name of the original file" << std::endl;
    std::cout << "   -i <input_file>            input_file must be the file the index was built for" << std::endl;
    std::cout << "                              (required for locate_rlzsa_bin_search and the -c option)" << std::endl;
    std::cout << "   -c                         checks correctness of each pattern occurrence on <input_file>" << std::endl;
    std::cout << "   -o <output_file>           write pattern occurrences to this file (ASCII)" << std::endl;
    std::cout << "   -t <arg>                   maximum number of occurrences to report per pattern. (Def. all)" << std::endl;
    std::cout << "   <index_file>               index file (with extension .move-r)" << std::endl;
    //std::cout << "   <patterns_file>            file in pizza&chili format containing the patterns" << std::endl;
    std::cout << "   <patterns_file>            file in FASTA format containing the patterns" << std::endl;
    exit(0);
}

void parse_args(char** argv, int argc, int& ptr)
{
    std::string s = argv[ptr];
    ptr++;

    if (s == "-c") {
        check_correctness = true;
    } else if (s == "-m") {
        if (ptr >= argc - 1) help("error: missing parameter after -o option.");
        std::string path_m_file = argv[ptr++];
        mf.open(path_m_file, std::filesystem::exists(path_m_file) ? std::ios::app : std::ios::out);
        if (!mf.good()) help("error: cannot open nor create <m_file>");
        name_text_file = argv[ptr++];
    } else if (s == "-i") {
        if (ptr >= argc - 1) help("error: missing parameter after -i option.");
        path_input_file = argv[ptr++];
        input_file.open(path_input_file);
        if (!input_file.good()) help("error: cannot open <input_file>");
    } else if (s == "-o") {
        if (ptr >= argc - 1) help("error: missing parameter after -o option.");
        output_occurrences = true;
        path_output_file = argv[ptr++];
    } else if (s == "-t") {
        if (ptr >= argc - 1) help("error: missing parameter after -t option.");
        max_occ = std::stoull(argv[ptr++]);
    } else {
        help("error: unrecognized '" + s + "' option");
    }
}

template <typename pos_t, move_r_support support>
void measure_locate()
{
    std::cout << "loading the index" << std::endl;
    //auto time = now();
    using idx_t = move_r<support, char, pos_t>;
    idx_t index;
    index.load(index_file);
    //time = log_runtime(time);
    index_file.close();
    std::cout << std::endl;
    //index.log_data_structure_sizes();
    /*
    if (support == _locate_rlzsa_bin_search || check_correctness)
    {
        std::cerr << "it sould not be entering here!" << std::endl;
        exit(1);
    }*/

    // --- RESET MEMORY COUNTER ---
    // compute index size + basic program overhead.
    size_t index_memory_bytes = malloc_count_current();

    // --- PHASE 1: LOAD DATA ---
    // read patterns into memory first to avoid disk I/O latency.
    std::ofstream   output  (path_patterns_file+".mover_occs");
    std::ofstream   logfile (path_patterns_file+".mover_locate_log");

    struct QueryData {
        std::string header;
        std::string pattern;
    };
    std::vector<QueryData> queries;

    using RangeType = decltype(index.count_and_get_range(""));
    using Step1ReturnType = std::pair<RangeType,int64_t>;
    std::vector<Step1ReturnType> step1_results;

    std::string line, header;
    uint64_t i = 0, tot_char = 0;

    while(std::getline(patterns_file, line)) 
    {
        if(i % 2 != 0) {
            queries.push_back({header, line});
            tot_char += line.size();
        } else {
            header = line;
        }
        i++;
    }
    patterns_file.close();

    // reserve space to prevent reallocation during timing
    step1_results.reserve(queries.size());

    uint64_t tot_pattern = queries.size();
    uint64_t tot_soccs = 0;
    double tot_primary_duration = 0;
    double tot_secondary_duration = 0;

    // --- PHASE 2: RUN STEP 1 (Primary occurrences) ---
    for(const auto& q : queries)
    {
        auto start = std::chrono::high_resolution_clock::now();
        auto res = index.count_and_get_range(q.pattern);
        int64_t pocc = std::get<1>(res) >= std::get<0>(res) ? 
                       index.primary_occ(std::get<2>(res),std::get<3>(res)) : -1;
        std::chrono::duration<double> duration = 
                std::chrono::high_resolution_clock::now() - start;

        tot_primary_duration += duration.count();
        step1_results.push_back(std::make_pair(res,pocc));
    }

    // --- RESET MEMORY COUNTER ---
    size_t memory_after_step_1 = malloc_count_current();
    malloc_count_reset_peak();

    // --- PHASE 3: RUN STEP 2 (Secondary occurrences) ---
    // Step 2 enumerates the secondary occurences.
    for(size_t k = 0; k < queries.size(); ++k)
    {
        auto& res = step1_results[k]; 
        int64_t pocc = res.second;

        auto start = std::chrono::high_resolution_clock::now();
        auto soccs = index.locate_secondary_occs(res.first, max_occ-1);
        std::chrono::duration<double> duration = 
                std::chrono::high_resolution_clock::now() - start;

        tot_secondary_duration += duration.count();
        tot_soccs += soccs.size();

        // write output outside the timer
        output << queries[k].header << std::endl;
        if( pocc != -1 )
        {   
            output << pocc << " ";
            for(auto& occ : soccs)
                output << occ << " ";
            output << std::endl;
        }
        else { output << std::endl; }
    }

    output.close();

    // --- PHASE 4: LOG DATA ---
    std::cout << "No. occurences found = "      << tot_soccs + tot_pattern              << "\n"
              << "No. processed patterns = "    << tot_pattern                          << "\n"
              << "No. processed characters = "  << tot_char                             << "\n" 
              << "Index size (bytes) = "        << index_memory_bytes                   << "\n"
              << "Memory peak secondary occ.s (bytes) = "     << malloc_count_peak() -
                                          (memory_after_step_1 - index_memory_bytes)    << "\n"
              << "Time to find the primary occ. (sec) = "     << tot_primary_duration   << "\n"
              << "Time to find the secondary occ.s (sec) = "  << tot_secondary_duration << "\n"
              << "Time per primary occurrence (ns) = "        << 
                                         (tot_primary_duration/tot_pattern) * 1000000000 << "\n"
              << "Time per secondary occurrence (ns) = "      << 
                                         (tot_secondary_duration/tot_soccs) * 1000000000 << std::endl;

    logfile   << "No. occurences found = "      << tot_soccs + tot_pattern              << "\n"
              << "No. processed patterns = "    << tot_pattern                          << "\n"
              << "No. processed characters = "  << tot_char                             << "\n" 
              << "Index size (bytes) = "        << index_memory_bytes                   << "\n"
              << "Memory peak secondary occ.s (bytes) = "     << malloc_count_peak() -
                                          (memory_after_step_1 - index_memory_bytes)    << "\n"
              << "Time to find the primary occ. (sec) = "     << tot_primary_duration   << "\n"
              << "Time to find the secondary occ.s (sec) = "  << tot_secondary_duration << "\n"
              << "Time per primary occurrence (ns) = "        << 
                                         (tot_primary_duration/tot_pattern) * 1000000000 << "\n"
              << "Time per secondary occurrence (ns) = "      << 
                                         (tot_secondary_duration/tot_soccs) * 1000000000 << std::endl;

    logfile.close();
}

template <typename pos_t, move_r_support support>
void measure_primary()
{
    std::cout << "loading the index" << std::endl;
    //auto time = now();
    using idx_t = move_r<support, char, pos_t>;
    idx_t index;
    index.load(index_file);
    //time = log_runtime(time);
    index_file.close();
    std::cout << std::endl;

    // --- RESET MEMORY COUNTER ---
    // compute index size + basic program overhead.
    size_t index_memory_bytes = malloc_count_current();

    // --- PHASE 1: LOAD DATA ---
    std::ofstream   logfile (path_patterns_file+".mover_locate_log");

    struct QueryData {
        std::string header;
        std::string pattern;
    };
    std::vector<QueryData> queries;

    std::string line, header;
    uint64_t i = 0, tot_char = 0;

    while(std::getline(patterns_file, line)) 
    {
        if(i % 2 != 0) {
            queries.push_back({header, line});
            tot_char += line.size();
        } else {
            header = line;
        }
        i++;
    }
    patterns_file.close();

    uint64_t tot_pattern = queries.size();
    double  tot_primary_duration = 0;
    uint64_t tot_primary_found = 0;

    // --- PHASE 2: RUN STEP 1 (Primary occurrences) ---
    for(const auto& q : queries)
    {
        auto start = std::chrono::high_resolution_clock::now();
        auto res = index.count_and_get_occ(q.pattern);
        int64_t pocc = res.first.second >= res.first.first ? 
                       res.second : -1;
        std::chrono::duration<double> duration = 
                std::chrono::high_resolution_clock::now() - start;

        tot_primary_duration += duration.count();
        tot_primary_found += (pocc == -1 ? 0 : 1);
    }

    // --- PHASE 4: LOG DATA ---
    std::cout << "No. primary occs found = "    << tot_primary_found                    << "\n"
              << "No. processed patterns = "    << tot_pattern                          << "\n"
              << "No. processed characters = "  << tot_char                             << "\n" 
              << "Index size (bytes) = "        << index_memory_bytes                   << "\n"
              << "Time to find the primary occ. (sec) = "     << tot_primary_duration   << "\n"
              << "Time per primary occurrence (ns) = " << std::fixed << std::setprecision(2) << 
                                        (tot_primary_duration/tot_pattern) * 1000000000 << std::endl;

    logfile   << "No. primary occs found = "    << tot_primary_found                    << "\n"
              << "No. processed patterns = "    << tot_pattern                          << "\n"
              << "No. processed characters = "  << tot_char                             << "\n" 
              << "Index size (bytes) = "        << index_memory_bytes                   << "\n"
              << "Time to find the primary occ. (sec) = "     << tot_primary_duration   << "\n"
              << "Time per primary occurrence (ns) = " << std::fixed << std::setprecision(2) << 
                                        (tot_primary_duration/tot_pattern) * 1000000000 << std::endl;

    logfile.close();
}

template <typename pos_t, move_r_support support>
void measure_secondary()
{
    std::cout << "loading the index" << std::endl;
    //auto time = now();
    using idx_t = move_r<support, char, pos_t>;
    idx_t index;
    index.load(index_file);
    index_file.close();

    // --- RESET MEMORY COUNTER ---
    // compute index size + basic program overhead.
    size_t index_memory_bytes = malloc_count_current();

    // --- PHASE 1: LOAD DATA ---
    // read patterns into memory first to avoid disk I/O latency.
    std::ofstream   logfile (path_patterns_file+".mover_locate_log");

    struct QueryData {
        std::string header;
        std::string pattern;
    };
    std::vector<QueryData> queries;

    std::string line, header;
    uint64_t i = 0, tot_char = 0;

    while(std::getline(patterns_file, line)) 
    {
        if(i % 2 != 0) {
            queries.push_back({header, line});
            tot_char += line.size();
        } else {
            header = line;
        }
        i++;
    }
    patterns_file.close();

    malloc_count_reset_peak();

    uint64_t tot_pattern = queries.size();
    uint64_t tot_soccs = 0;
    double tot_secondary_duration = 0;

    // --- PHASE 2: RUN QUERIES (Secondary occurrences) ---
    for(const auto& q : queries)
    {
        auto start = std::chrono::high_resolution_clock::now();

        auto res = index.count_and_get_range(q.pattern);
        int64_t pocc = std::get<1>(res) >= std::get<0>(res) ? 
                       index.primary_occ(std::get<2>(res),std::get<3>(res)) : -1;

        std::vector<pos_t> soccs{};
        if(pocc != -1)
            soccs = index.locate_secondary_occs(res, max_occ-1);

        std::chrono::duration<double> duration = 
                std::chrono::high_resolution_clock::now() - start;

        tot_secondary_duration += duration.count();
        tot_soccs += soccs.size() + (pocc == -1 ? 0 : 1);
    }

    // --- PHASE 4: LOG DATA ---
    std::cout << "No. occurences found = "      << tot_soccs                            << "\n"
              << "No. processed patterns = "    << tot_pattern                          << "\n"
              << "No. processed characters = "  << tot_char                             << "\n" 
              << "Index size (bytes) = "        << index_memory_bytes                   << "\n"
              << "Memory peak secondary occ.s (bytes) = "     << malloc_count_peak()    << "\n"
              << "Time to find the secondary occ.s (sec) = "  << tot_secondary_duration << "\n"
              << "Time per secondary occurrence (ns) = "      << 
                                         (tot_secondary_duration/tot_soccs) * 1000000000 << std::endl;

    logfile   << "No. occurences found = "      << tot_soccs + tot_pattern              << "\n"
              << "No. processed patterns = "    << tot_pattern                          << "\n"
              << "No. processed characters = "  << tot_char                             << "\n" 
              << "Index size (bytes) = "        << index_memory_bytes                   << "\n"
              << "Memory peak secondary occ.s (bytes) = "     << malloc_count_peak()    << "\n"
              << "Time to find the secondary occ.s (sec) = "  << tot_secondary_duration << "\n"
              << "Time per secondary occurrence (ns) = "      << 
                                         (tot_secondary_duration/tot_soccs) * 1000000000 << std::endl;

    logfile.close();
}

int main(int argc, char** argv)
{
    if (argc < 3) help("");
    while (ptr < argc - 2) parse_args(argv, argc, ptr);

    path_index_file = argv[ptr];
    path_patterns_file = argv[ptr + 1];

    index_file.open(path_index_file);
    patterns_file.open(path_patterns_file);

    if (!index_file.good()) help("error: could not read <index_file>");
    if (!patterns_file.good()) help("error: could not read <patterns_file>");

    if (output_occurrences) {
        output_file.open(path_output_file);
        if (!output_file.good()) help("error: could not create <output_file>"); 
    }

    bool is_64_bit;
    index_file.read((char*) &is_64_bit, 1);
    move_r_support _support;
    index_file.read((char*) &_support, sizeof(move_r_support));
    index_file.seekg(0, std::ios::beg);

    if(is_64_bit)
        std::cout << "running in 64 bit environment." << std::endl;
    else
        std::cout << "running in 32 bit environment." << std::endl;

    if (_support == _count) {
        std::cout << "error: this index does not support locate" << std::endl;
        exit(0);
    }
    else if (_support == _locate_one) {
        if (is_64_bit) {
            measure_primary<uint64_t, _locate_one>();
        } else {
            measure_primary<uint32_t, _locate_one>();
        }
    }
    else if (_support == _locate_move) {
        if (is_64_bit) {
            measure_secondary<uint64_t, _locate_move>();
        } else {
            measure_secondary<uint32_t, _locate_move>();
        }
    }
    else {
        std::cout << "error: we are only testing the move and locate_one variant!" << std::endl;
        exit(1);
    }
    /*
    else if (_support == _locate_rlzsa) {
        if (is_64_bit) {
            measure_locate<uint64_t, _locate_rlzsa>();
        } else {
            measure_locate<uint32_t, _locate_rlzsa>();
        }
    } else if (_support == _locate_rlzsa_bin_search) {
        if (is_64_bit) {
            measure_locate<uint64_t, _locate_rlzsa_bin_search>();
        } else {
            measure_locate<uint32_t, _locate_rlzsa_bin_search>();
        }
    } else if (_support == _locate_lzendsa) {
        if (is_64_bit) {
            measure_locate<uint64_t, _locate_lzendsa>();
        } else {
            measure_locate<uint32_t, _locate_lzendsa>();
        }
    }
    */

    patterns_file.close();
    if (output_occurrences)
        output_file.close();
}