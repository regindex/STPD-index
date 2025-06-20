set(RLZ_SOURCES RLZ_DNA.hpp RLZ_DNA_sux.hpp)

add_library(RLZ OBJECT ${RLZ_SOURCES})
target_link_libraries(RLZ PUBLIC common elias_fano sdsl divsufsort divsufsort64)