include(FetchContent)

######### malloc_count #########

  FetchContent_Declare(
    malloc_count
    GIT_REPOSITORY https://github.com/bingmann/malloc_count.git
  )

  FetchContent_GetProperties(malloc_count)
  if(NOT malloc_count_POPULATED)
    FetchContent_Populate(malloc_count)
  endif()

  add_library(malloc_count STATIC
    ${malloc_count_SOURCE_DIR}/malloc_count.c
    ${malloc_count_SOURCE_DIR}/malloc_count.h
  )
  target_include_directories(malloc_count PUBLIC "${malloc_count_SOURCE_DIR}")
  target_link_libraries(malloc_count PUBLIC dl)

  add_library(memprofile OBJECT
    ${malloc_count_SOURCE_DIR}/memprofile.h
  )
  target_include_directories(memprofile PUBLIC "${malloc_count_SOURCE_DIR}")

######### sdsl-lite #########

  FetchContent_Declare(
    sdsl
    GIT_REPOSITORY https://github.com/simongog/sdsl-lite
    GIT_TAG v2.1.1
  )

  set(GENERATE_DOC OFF CACHE BOOL "Disable Doxygen generation for SDSL" FORCE)
  set(SDSL_BUILD_EXAMPLES OFF CACHE BOOL "Do not build SDSL examples" FORCE)
  set(SDSL_BUILD_TESTS OFF CACHE BOOL "Do not build SDSL tests" FORCE)

  FetchContent_MakeAvailable(sdsl)

######### divsufsort #########

  #FetchContent_Declare(
  #  divsufsort
  #  GIT_REPOSITORY https://github.com/simongog/libdivsufsort.git
  #  GIT_TAG        2.0.1
  #)

  #set(BUILD_SHARED_LIBS OFF CACHE BOOL "Do not build a shared library for libdivsufsort")
  #set(BUILD_EXAMPLES OFF CACHE BOOL "Do not build libdivsufsort example")
  #set(BUILD_DIVSUFSORT64 ON CACHE BOOL "Build libdivsufsort in 64-bits mode")

  #FetchContent_MakeAvailable(divsufsort)

######### gsacak #########

  FetchContent_Declare(
    gsacak
    GIT_REPOSITORY https://github.com/felipelouza/gsa-is.git
  )

  FetchContent_GetProperties(gsacak)
  if(NOT gsacak_POPULATED)
      FetchContent_Populate(gsacak)
  endif()

  add_library(gsacak OBJECT
    ${gsacak_SOURCE_DIR}/gsacak.c
    ${gsacak_SOURCE_DIR}/gsacak.h
  )
  target_include_directories(gsacak PUBLIC "${gsacak_SOURCE_DIR}")

  add_library(gsacak64 OBJECT
    ${gsacak_SOURCE_DIR}/gsacak.c
    ${gsacak_SOURCE_DIR}/gsacak.h
  )
  target_include_directories(gsacak64 PUBLIC "${gsacak_SOURCE_DIR}")
  target_compile_options(gsacak64 PUBLIC -DM64)

######### sux #########

  add_library(sux OBJECT
    sux/Rank.hpp
    sux/SimpleSelectHalf.hpp
    sux/SimpleSelectZeroHalf.hpp
  )
  target_include_directories(sux PUBLIC sux)