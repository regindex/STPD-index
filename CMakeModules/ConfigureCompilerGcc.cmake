# ##############################################################################
# Gcc Compiler configuration
# ##############################################################################

# Set C compiler standard
set(CMAKE_C_STANDARD 11)    
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
# Set C++ compiler standard
set(CMAKE_CXX_STANDARD 20) 
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add the basic compiler options
add_compile_options("-Wall")
add_compile_options("-Wextra")
add_compile_options("-Wcomment")

# Compiler flags for all different configuration
set(DEBUG_FLAGS "-fsanitize=address -fno-omit-frame-pointer") # debug mode
set(RELEASE_FLAGS "-march=native -funroll-loops -fomit-frame-pointer -flto -w") # release mode
set(RELDEBINFO_FLAGS "-fno-omit-frame-pointer") # release wirh debug info
set(RELMINSIZ_FLAGS "-march=native -funroll-loops -fomit-frame-pointer -flto -w") # release with minimum binaries size

# Apply flags globally
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DEBUG_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${RELEASE_FLAGS}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${RELDEBINFO_FLAGS}")
set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} ${RELMINSIZ_FLAGS}")

# Print configuration info
message(STATUS "C standard: ${CMAKE_C_STANDARD}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Debug flags: ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "Release flags: ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "RelWithDebInfo flags: ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
message(STATUS "MinSizeRel flags: ${CMAKE_CXX_FLAGS_MINSIZEREL}")
