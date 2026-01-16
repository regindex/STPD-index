### Cmake initialization ###

  cmake_minimum_required(VERSION "3.14.0")
  set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")

  cmake_policy(SET CMP0048 NEW)
  set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)

### About this project ###

  project(STPD-index LANGUAGES C CXX)
  SET(VERSION_MAJOR "2")
  SET(VERSION_MINOR "0")
  SET(VERSION_PATCH "0")
  SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

### Set the build type ###

  if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  endif()

  message("building in ${CMAKE_BUILD_TYPE} mode")

### Set environment and compiler ###

  set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install)

  find_package(Git REQUIRED)
  if(GIT_FOUND)
      message("git found: ${GIT_EXECUTABLE}")
  else()
      message(WARNING "git not found. Cloning of submodules will not work.")
  endif()

  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
      # using Clang
      message("compiler: configuring Clang compiler")
      include(ConfigureCompilerClang)
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
      # using GCC
      message("compiler: configuring gcc compiler")
      include(ConfigureCompilerGcc)
  else ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
      message(FATAL_ERROR "Only the compiler gcc and clang are supported")
  endif()

### CMake options ###

  option(STPD_BUILD_CLI "Build Command Line Interface (CLI) executables" ON)
  option(STPD_BUILD_BENCH_CLI "Build Benchmark Executables" ON)

### Configure external dependencies ###

  add_subdirectory(external)

### Configure internal resources ###

  set(CMAKE_INSTALL_INCLUDEDIR "include")
  add_subdirectory(include)

  add_library(stpd INTERFACE)
  target_include_directories(stpd INTERFACE "${CMAKE_SOURCE_DIR}/include/")
  target_link_libraries(stpd INTERFACE
    sdsl malloc_count sux divsufsort divsufsort64 common bitvector elias_fano
    phi move_phi text_oracle binary_search text_index
  )

### Configure CLI executables ###

  if(STPD_BUILD_CLI)
    add_subdirectory(CLI)
  endif()

### Configure BENCHMARK FILES executables ###

  if(STPD_BUILD_BENCH_CLI)
    add_subdirectory(benchmark)
  endif()