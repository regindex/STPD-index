set(EF_SOURCES elias_fano_sux.hpp)

add_library(elias_fano OBJECT ${EF_SOURCES})
target_link_libraries(elias_fano sux)