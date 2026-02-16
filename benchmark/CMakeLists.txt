include(FetchContent)

message(STATUS "Configuring competitor software included in the benchmarks...")

######### r-index #########

  #FetchContent_Declare(
  #  r_index
  #  GIT_REPOSITORY https://github.com/nicolaprezza/r-index.git
  #  GIT_TAG master
  #)

  #FetchContent_GetProperties(r_index)
  #if(NOT r_index_POPULATED)
  #  FetchContent_Populate(r_index)
  #endif()

  #message(STATUS "r-index downloaded to: ${r_index_SOURCE_DIR}")

  #file(COPY ${CMAKE_SOURCE_DIR}/benchmark/r-index_patched/ri-locate-benchmark.cpp
  # DESTINATION ${r_index_SOURCE_DIR})

  #file(COPY ${CMAKE_SOURCE_DIR}/benchmark/r-index_patched/r_index.hpp
  #   DESTINATION ${r_index_SOURCE_DIR}/internal)

  message(STATUS "Building the r-index...")

  add_executable(ri_locate
    ${CMAKE_SOURCE_DIR}/benchmark/r-index/ri-locate-benchmark.cpp
  )

  target_include_directories(ri_locate PUBLIC
      ${CMAKE_SOURCE_DIR}/benchmark/r-index
      ${CMAKE_SOURCE_DIR}/benchmark/r-index/internal
  )

  target_link_libraries(ri_locate PUBLIC malloc_count sdsl divsufsort divsufsort64)

  add_executable(ri_build
    ${CMAKE_SOURCE_DIR}/benchmark/r-index/ri-build.cpp
  )

  target_include_directories(ri_build PUBLIC
      ${CMAKE_SOURCE_DIR}/benchmark/r-index
      ${CMAKE_SOURCE_DIR}/benchmark/r-index/internal
  )

  target_link_libraries(ri_build PUBLIC malloc_count sdsl divsufsort divsufsort64)

######### Move-r #########

##-- Move-r dependencies --##

  find_package(OpenMP REQUIRED)

  add_subdirectory("${CMAKE_SOURCE_DIR}/benchmark/Move-r/external/ips4o")
  set_target_properties(ips4o_example PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)

  add_subdirectory("${CMAKE_SOURCE_DIR}/benchmark/Move-r/external/concurrentqueue/")

  set(LIBSAIS_USE_OPENMP ON CACHE BOOL "" FORCE)
  add_subdirectory("${CMAKE_SOURCE_DIR}/benchmark/Move-r/external/libsais/")

  add_library(malloc_count_ INTERFACE)
  target_sources(malloc_count_ INTERFACE "${CMAKE_SOURCE_DIR}/benchmark/Move-r/external/malloc_count/malloc_count.c")
  target_include_directories(malloc_count_ INTERFACE "${CMAKE_SOURCE_DIR}/benchmark/Move-r/external/malloc_count/")
  target_link_libraries(malloc_count_ INTERFACE dl)

  add_subdirectory("${CMAKE_SOURCE_DIR}/benchmark/Move-r/external/gtl/")

  add_subdirectory("${CMAKE_SOURCE_DIR}/benchmark/Move-r/external/sparse-map/")

  add_library(emhash INTERFACE)
  target_include_directories(emhash INTERFACE "${CMAKE_SOURCE_DIR}/benchmark/Move-r/external/emhash/")

  add_library(sux_ INTERFACE)
  target_include_directories(sux_ INTERFACE "${CMAKE_SOURCE_DIR}/benchmark/Move-r/external/sux/")

  add_subdirectory("${CMAKE_SOURCE_DIR}/benchmark/Move-r/external/ordered/")

  add_library(rmq INTERFACE)
  target_include_directories(rmq INTERFACE "${CMAKE_SOURCE_DIR}/benchmark/Move-r/external/rmq/include/")

  add_subdirectory("${CMAKE_SOURCE_DIR}/benchmark/Move-r/external/unordered_dense/")

  ##-- setup Move-r library and executables --##

  add_library(move_r INTERFACE) 
  target_include_directories(move_r INTERFACE "${CMAKE_SOURCE_DIR}/benchmark/Move-r/include/")
  target_link_libraries(move_r INTERFACE
    OpenMP::OpenMP_CXX sdsl libsais ips4o gtl sparse_map malloc_count_
    concurrentqueue emhash sux_ ordered rmq divsufsort divsufsort64
  )

  add_executable(move-r-build  ${CMAKE_SOURCE_DIR}/benchmark/Move-r/cli/move-r/move-r-build.cpp)
  add_executable(move-r-locate ${CMAKE_SOURCE_DIR}/benchmark/Move-r/cli/move-r/move-r-locate.cpp)

  target_link_libraries(move-r-build  PRIVATE move_r)
  target_link_libraries(move-r-locate PRIVATE move_r)

######### PA-index #########

  add_executable(pa_build
    pa-index/build_store_pa_index.cpp
  )

  target_link_libraries(pa_build PUBLIC malloc_count text_oracle)

  add_executable(pa_locate
    pa-index/locate.cpp
  )

  target_link_libraries(pa_locate PUBLIC malloc_count text_oracle)

######### PT-index #########

  add_executable(pt_build
    prefix-tree-index/build_pt_index.cpp
  )

  target_link_libraries(pt_build PUBLIC malloc_count text_oracle)

  add_executable(pt_locate
    prefix-tree-index/locate.cpp
  )

  target_link_libraries(pt_locate PUBLIC malloc_count text_oracle)

######### Suffixient-array #########

  add_executable(build_suffixient
    suffixient-array/compute_suffixient_set.cpp
  )

  target_link_libraries(build_suffixient PRIVATE stpd)

  add_executable(sA-build
    suffixient-array/build_sA-index.cpp
  )

  target_link_libraries(sA-build PRIVATE stpd)

  add_executable(sA_locate
    suffixient-array/locate.cpp
  )

  target_link_libraries(sA_locate PRIVATE stpd)

######### script to build the indexes variants #########

  #add_executable(benchmark_build
  #  build_index_variants.cpp
  #)

  #target_link_libraries(benchmark_build PRIVATE stpd)

  add_executable(sampling_build
    ${CMAKE_SOURCE_DIR}/CLI/path-decomp-src/STPD_small_files.cpp
  )

  target_link_libraries(sampling_build PRIVATE stpd)

######### benchmark script #########

  #file(COPY ${CMAKE_SOURCE_DIR}/benchmark/testing_script.py
  #   DESTINATION ${CMAKE_CURRENT_BINARY_DIR})