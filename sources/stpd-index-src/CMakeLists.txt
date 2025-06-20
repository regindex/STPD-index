add_executable(build_store_stpd_index build_store_stpd_index.cpp)
#target_link_libraries(build_store_stpd_index bitvectors common RLZ phi_functions malloc_count sdsl divsufsort divsufsort64 pthread) 
target_link_libraries(build_store_stpd_index PUBLIC RLZ stpd_array phi_functions malloc_count) 

add_executable(locate locate.cpp)
#target_link_libraries(locate bitvectors common RLZ phi_functions malloc_count sdsl divsufsort divsufsort64 pthread) 
target_link_libraries(locate PUBLIC RLZ stpd_array phi_functions malloc_count) 

add_executable(locate64 locate.cpp)
target_link_libraries(locate64 PUBLIC RLZ stpd_array phi_functions malloc_count)
target_compile_options(locate64 PUBLIC "-DM64")

#add_executable(test test.cpp)
#target_link_libraries(test elias_fano common) 