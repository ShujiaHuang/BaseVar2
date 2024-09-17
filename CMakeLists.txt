# how to write CMakeList.txt: https://subingwen.cn/cmake/CMake-prime
cmake_minimum_required (VERSION 3.12.0)

# Project information
project (BaseVar VERSION 1.1.0)

include(ExternalProject)

# Set the C++ standard required for the project
set(CMAKE_CXX_STANDARD 17)  # C++ 17
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3 -fPIC")  # default parameters for building.
message(STATUS "CXX flags are: ${CMAKE_CXX_FLAGS}")

if (DEFINED CFLAGS)
    set(MYCFLAGS "${CFLAGS}")
    message (STATUS "Custom CFLAGS: ${MYCFLAGS}")
endif()

if (DEFINED LDFLAGS)
  set(MYLDFLAGS "${LDFLAGS}")
  message (STATUS "Custom LDFLAGS: ${MYLDFLAGS}")
endif()

set (STATIC_DIR "" CACHE STRING "Build in 'static' mode and include libraries in this directory.")

#############################################################################
## Compile htslib
#############################################################################
if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/htslib/Makefile)
    message(FATAL_ERROR "htslib not found. Run 'git submodule update --init' to fetch all submodules.")
endif()

# add_custom_target(
ExternalProject_Add(
    project_htslib
    BUILD_IN_SOURCE 1
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/htslib
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/htslib
    CONFIGURE_COMMAND cp -a ${CMAKE_CURRENT_SOURCE_DIR}/htslib ${CMAKE_CURRENT_BINARY_DIR}/ COMMAND autoreconf -i COMMAND ${CMAKE_CURRENT_BINARY_DIR}/htslib/configure --disable-libcurl --disable-gcs
        "CFLAGS=${MYCFLAGS} -g -Wall -O3 ${CMAKE_C_FLAGS}"
        "LDFLAGS=${MYLDFLAGS}"
        "CC=${CMAKE_C_COMPILER}"
    BUILD_COMMAND $(MAKE) -C ${CMAKE_CURRENT_BINARY_DIR}/htslib
)

add_library(htslib STATIC IMPORTED)

include_directories(${PROJECT_SOURCE_DIR}/htslib)      # include the header of htslib
link_libraries(${PROJECT_SOURCE_DIR}/htslib/libhts.a)  # link htslib static lib: libhts.a
# set(htslib_location ${CMAKE_CURRENT_BINARY_DIR}/htslib/libhts.a)
# message(STATUS "htslib target location is ${htslib_location}")
# set_property(TARGET htslib PROPERTY IMPORTED_LOCATION ${htslib_location})
# add_dependencies(htslib project_htslib)

# LOCAL binaries have static GCC, PREBUILT are all static
if (STATIC_DIR STREQUAL "")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
else ()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")
endif ()

#############################################################################
##  DEPENDENCY END
#############################################################################
message(STATUS "Libraries: ${COMMON_SYS_LIB}")
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

# Find all relevent source codes
file(GLOB SRC_LIST ${PROJECT_SOURCE_DIR}/src/*.cpp)
set(COMMON_SYS_LIB pthread z bz2 m lzma curl)

# Create executable basevar
add_executable(basevar ${SRC_LIST})
if (APPLE)
    target_link_libraries(basevar ${COMMON_SYS_LIB})
else ()
    target_link_libraries(basevar ${COMMON_SYS_LIB} ssl crypto)
endif()
 