set(STPD_SOURCES stpd_array_binary_search.hpp stpd_array_binary_search_opt.hpp stpd_array_binary_search_opt_v2.hpp)

add_library(stpd_array OBJECT ${STPD_SOURCES})
target_link_libraries(common text_oracles bitvectors elias_fano sdsl divsufsort divsufsort64)