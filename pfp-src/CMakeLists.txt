add_executable(pscan pscan.cpp utils.c xerrors.c)
target_link_libraries(pscan pthread) 

add_executable(bwtparse bwtparse.c utils.c)
target_link_libraries(bwtparse gsacak) 

add_executable(pfbwt pfbwt.cpp utils.c xerrors.c)
target_link_libraries(pfbwt pthread gsacak) 