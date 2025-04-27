# ##############################################################################
# Clang compiler configuration
# ##############################################################################

set(CMAKE_C_STANDARD 11)    
set(CMAKE_CXX_STANDARD 14)   

# Add the basic compiler options
# add_compile_options("-Werror")
add_compile_options("-Wall")
add_compile_options("-Wextra")
add_compile_options("-Wcomment")


# Add the basic compiler options for debug version
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb3 -fsanitize=address -fno-omit-frame-pointer")
# Add the basic compiler options for release version
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ansi -march=native -funroll-loops -O3 -DNDEBUG -w")
