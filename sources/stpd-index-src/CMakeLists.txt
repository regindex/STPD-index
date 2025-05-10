add_executable(build_store_stpd_index build_store_stpd_index.cpp)
target_link_libraries(build_store_stpd_index bitvectors common RLZ phi_functions malloc_count sdsl divsufsort divsufsort64 pthread) 

add_executable(locate locate.cpp)
target_link_libraries(locate bitvectors common RLZ phi_functions malloc_count sdsl divsufsort divsufsort64 pthread) 