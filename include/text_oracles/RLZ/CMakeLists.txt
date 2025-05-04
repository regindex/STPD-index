set(RLZ_SOURCES RLZ_DNA.hpp)

add_library(RLZ OBJECT ${RLZ_SOURCES})
target_link_libraries(common sdsl divsufsort divsufsort64)