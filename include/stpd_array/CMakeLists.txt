set(STPD_SOURCES stpd_array_binary_search.hpp stpd_array_binary_search_opt.hpp)

add_library(stpd_array OBJECT ${STPD_SOURCES})
target_link_libraries(stpd_array PUBLIC common text_oracles bitvectors)