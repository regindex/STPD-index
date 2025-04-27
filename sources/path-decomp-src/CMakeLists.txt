add_executable(stpd_small STPD_small_files.cpp)
target_link_libraries(stpd_small sdsl divsufsort divsufsort64) 

#add_executable(colex+ compute_stpd_sampling.cpp)
#target_link_libraries(pscan pthread) 

#add_executable(bwtparse bwtparse.c utils.c)
#target_link_libraries(bwtparse gsacak) 

#add_executable(pfbwt pfbwt.cpp utils.c xerrors.c)
#target_link_libraries(pfbwt pthread gsacak) 