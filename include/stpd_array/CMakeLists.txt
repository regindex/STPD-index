set(STPD_SOURCES stpd_array_binary_search.hpp)

add_library(stpd_array OBJECT ${STPD_SOURCES})
target_link_libraries(common text_oracles bitvectors sdsl divsufsort divsufsort64)