# ##############################################################################
# Clang compiler configuration
# ##############################################################################

set(CMAKE_C_STANDARD 11)    
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 14)  
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the basic compiler options
# add_compile_options("-Werror")
add_compile_options("-Wall")
add_compile_options("-Wextra")
add_compile_options("-Wcomment")


# Add the basic compiler options for debug version
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb3 -O0 -DDEBUG")
# Add the basic compiler options for release version
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ansi -march=native -funroll-loops -O3 -DNDEBUG -w")

message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
