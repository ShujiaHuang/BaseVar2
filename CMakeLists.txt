# how to write CMakeList.txt: https://subingwen.cn/cmake/CMake-primer
cmake_minimum_required (VERSION 3.12.0)

# Set the C++ standard required for the project
set(CMAKE_CXX_STANDARD 17)  # C++ 17
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set the project name
project (BASEVAR VERSION 2.2.0 LANGUAGES CXX)

# Add project metadata
set(PROJECT_AUTHOR "Shujia Huang")
set(PROJECT_EMAIL "hshujia@qq.com")
set(PROJECT_DESCRIPTION "BaseVar: Call variants and estimate allele frequency from ultra low-pass (<1x) WGS data")

# Configure version header
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/src/version.h.in
    ${CMAKE_CURRENT_SOURCE_DIR}/src/version.h
)

# Build options
option(STATIC_BUILD "Build statically linked executable" OFF)

# Set compiler flags
# -O3: optimization level 3
# -fPIC: generate position-independent code
# -Wl,-no_compact_unwind: disable compact unwind information for MacOS only
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fPIC")
if (APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,-no_compact_unwind")
endif()

# -------------------------------------------------------
# Static linking configuration
# -------------------------------------------------------
# Linux  : Fully static executable (-static), works with musl (Alpine) or glibc.
#          On glibc hosts, some features (e.g. DNS via NSS) may not work in
#          fully-static mode.  Building inside an Alpine (musl) container is
#          recommended for maximum portability.
# macOS  : Fully static executables are NOT supported by Apple.  We statically
#          link the C++ runtime (-static-libgcc -static-libstdc++) and all
#          third-party .a archives.  Only system frameworks (libSystem, etc.)
#          remain dynamically linked.
# -------------------------------------------------------
if(STATIC_BUILD)
    if(APPLE)
        # Apple Clang does NOT support -static-libgcc / -static-libstdc++.
        # Strategy: pass all third-party archives explicitly as .a files and
        # let the linker prefer static over dynamic for those.  Only Apple
        # system frameworks (/usr/lib/libSystem.B.dylib, etc.) will remain
        # dynamic — this is the best achievable portability on macOS.
        #
        # We intentionally leave CMAKE_EXE_LINKER_FLAGS unchanged here;
        # the static .a paths are handled below in target_link_libraries
        # by replacing the -l<name> flags with explicit .a paths.
        message(STATUS "STATIC_BUILD is ON (macOS) — linking third-party libs statically via explicit .a paths")
    else()
        # Linux static build: -static linker flag is set in the
        # target_link_libraries section below, after pkg-config detection.
        message(STATUS "STATIC_BUILD is ON (Linux) — producing a fully static basevar")
    endif()
endif()

# Ensure htslib is built before the main executable.
# For static builds on Linux, disable libcurl in htslib: BaseVar only reads
# local BAM/CRAM/FASTA files, so remote-URL support is not needed.
# This removes curl (and its transitive deps: brotli, zstd, psl, etc.)
# from the link graph, making the static build much simpler.
if(STATIC_BUILD AND NOT APPLE)
    set(HTSLIB_CONFIGURE_FLAGS --disable-libcurl)
else()
    set(HTSLIB_CONFIGURE_FLAGS "")
endif()

add_custom_target(
    Buildhts ALL
    COMMAND autoreconf -i && ./configure ${HTSLIB_CONFIGURE_FLAGS} && make
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

# Link libraries
if(STATIC_BUILD AND APPLE)
    # ---- macOS best-effort static build ----
    # Apple does not support fully static executables.
    # Strategy: statically link zlib, bz2, lzma via Homebrew .a archives.
    # Homebrew's libcurl.a requires ngtcp2/libicucore (QUIC/HTTP3+IDN) which
    # are complex to resolve as .a on macOS — use system curl dylib instead.
    #
    # Homebrew on Apple Silicon: /opt/homebrew/opt/<pkg>
    # Homebrew on Intel Macs:    /usr/local/opt/<pkg>
    find_library(ZLIB_STATIC   libz.a
        HINTS /opt/homebrew/opt/zlib/lib /usr/local/opt/zlib/lib
        REQUIRED NO_DEFAULT_PATH)
    find_library(BZ2_STATIC    libbz2.a
        HINTS /opt/homebrew/opt/bzip2/lib /usr/local/opt/bzip2/lib
        REQUIRED NO_DEFAULT_PATH)
    find_library(LZMA_STATIC   liblzma.a
        HINTS /opt/homebrew/opt/xz/lib /usr/local/opt/xz/lib
        REQUIRED NO_DEFAULT_PATH)

    message(STATUS "Static libs: zlib=${ZLIB_STATIC} bz2=${BZ2_STATIC} lzma=${LZMA_STATIC}")
    message(STATUS "Dynamic (system): curl, ssl, crypto")

    target_link_libraries(basevar
        ${HTSLIB}
        curl            # system /usr/lib/libcurl.4.dylib
        ${LZMA_STATIC}
        ${BZ2_STATIC}
        ${ZLIB_STATIC}
        pthread
        m
    )
elseif(STATIC_BUILD AND NOT APPLE)
    # ---- Linux fully static build (Alpine musl recommended) ----
    # htslib is configured with --disable-libcurl (see Buildhts above),
    # so curl and all its transitive deps (brotli, zstd, psl, etc.) are
    # not needed. Only the core compression libs + ssl are required.
    find_library(ZLIB_STATIC   NAMES libz.a     REQUIRED)
    find_library(BZ2_STATIC    NAMES libbz2.a   REQUIRED)
    find_library(LZMA_STATIC   NAMES liblzma.a  REQUIRED)
    find_library(SSL_STATIC    NAMES libssl.a   REQUIRED)
    find_library(CRYPTO_STATIC NAMES libcrypto.a REQUIRED)

    message(STATUS "Linux static libs:")
    message(STATUS "  zlib=${ZLIB_STATIC} bz2=${BZ2_STATIC} lzma=${LZMA_STATIC}")
    message(STATUS "  ssl=${SSL_STATIC} crypto=${CRYPTO_STATIC}")

    # Force static linker flag
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")

    target_link_libraries(basevar
        ${HTSLIB}
        ${SSL_STATIC}
        ${CRYPTO_STATIC}
        ${LZMA_STATIC}
        ${BZ2_STATIC}
        ${ZLIB_STATIC}
        pthread
        m
        dl
    )
else()
    # ---- Standard dynamic build ----
    target_link_libraries(basevar ${HTSLIB} ${COMMON_SYS_LIB})
endif()
