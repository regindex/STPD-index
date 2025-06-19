set(EF_SOURCES elias_fano_sux.hpp elias_fano_intlv.hpp)

add_library(elias_fano OBJECT ${EF_SOURCES})
target_link_libraries(elias_fano sux)