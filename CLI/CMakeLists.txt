add_executable(stpd-sampling  path-decomp-src/STPD_small_files.cpp)
add_executable(stpd-build     stpd-index-src/build_index.cpp)
add_executable(stpd-locate    stpd-index-src/locate.cpp)

target_link_libraries(stpd-sampling  PRIVATE stpd)
target_link_libraries(stpd-build     PRIVATE stpd)
target_link_libraries(stpd-locate    PRIVATE stpd)

# Set binaries location
set_target_properties(
    stpd-sampling stpd-build stpd-locate
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/CLI/"
)