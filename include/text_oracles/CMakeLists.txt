set(ORACLES_SOURCES bitpacked_text_oracle.hpp)

add_library(text_oracles OBJECT ${ORACLES_SOURCES})
target_link_libraries(common stpd_array text_oracles bitvectors sdsl divsufsort divsufsort64)

add_subdirectory(RLZ)
target_include_directories(RLZ PUBLIC RLZ)