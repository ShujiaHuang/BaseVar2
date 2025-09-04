# how to write CMakeList.txt: https://subingwen.cn/cmake/CMake-primer
cmake_minimum_required (VERSION 3.12.0)

# Set the C++ standard required for the project
set(CMAKE_CXX_STANDARD 17)  # C++ 17
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set the project name
project (BASEVAR VERSION 2.1.0 LANGUAGES CXX)

# Add project metadata
set(PROJECT_AUTHOR "Shujia Huang")
set(PROJECT_EMAIL "hshujia@qq.com")
set(PROJECT_DESCRIPTION "BaseVar: Call variants and estimate allele frequency from ultra low-pass (<1x) WGS data")

# Configure version header
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/src/version.h.in
    ${CMAKE_CURRENT_SOURCE_DIR}/src/version.h
)

# Set compiler flags
# -O3: optimization level 3
# -fPIC: generate position-independent code
# -Wl,-no_compact_unwind: disable compact unwind information for MacOS only
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fPIC")
if (APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,-no_compact_unwind")
endif()

# Ensure htslib is built before the main executable
add_custom_target(
    Buildhts ALL
    COMMAND autoreconf -i && ./configure && make
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/htslib
)

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/htslib)  # include the header of htslib
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src) 

# Define common system libraries
set(COMMON_SYS_LIB pthread z bz2 m lzma curl)
if(NOT APPLE)
    list(APPEND COMMON_SYS_LIB ssl crypto)
endif()

# Define htslib library
set(HTSLIB ${CMAKE_CURRENT_SOURCE_DIR}/htslib/libhts.a)

# Find all source codes
file(GLOB SRC_LIST 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/io/*.cpp
)

# create executable basevar
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)
add_executable(basevar ${SRC_LIST})
add_dependencies(basevar Buildhts)  # Ensure htslib is built before basevar
target_link_libraries(basevar ${HTSLIB} ${COMMON_SYS_LIB})
