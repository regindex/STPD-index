set(BITVECTORS_SOURCES elias_fano_bitvector.hpp succinct_bitvector.hpp)

add_library(bitvectors OBJECT ${BITVECTORS_SOURCES})
target_link_libraries(bitvectors PUBLIC common sdsl divsufsort divsufsort64)