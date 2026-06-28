# how to write CMakeList.txt: https://subingwen.cn/cmake/CMake-primer
cmake_minimum_required (VERSION 3.12.0)

# Set the C++ standard required for the project
set(CMAKE_CXX_STANDARD 17)  # C++ 17
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set the project name
project (BASEVAR VERSION 2.5.2 LANGUAGES CXX)

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
# Linux  : "Portable" partial-static executable.  We build on a glibc
#          host (Ubuntu in CI) and statically link:
#            * libstdc++ and libgcc      (-static-libstdc++ -static-libgcc)
#            * htslib + zlib + bzip2 + lzma + openssl  (.a archives)
#          glibc itself remains dynamic (the only sane choice on glibc
#          systems).  The result runs unmodified on any reasonably
#          modern Linux distro (Ubuntu 20.04+, Debian 11+, RHEL 8+,
#          CentOS 8+, Fedora 32+, ...).
#
#          We deliberately abandoned the previous Alpine-musl fully-
#          static build in v2.2.2: the resulting binary segfaulted
#          under multi-threaded workloads on glibc hosts due to
#          musl/glibc runtime incompatibilities (mallocng and pthread
#          behave differently on a glibc kernel under load).  Building
#          against the same libc the user runs on eliminates an entire
#          class of issues at once.
#
# macOS  : Fully static executables are NOT supported by Apple.  We
#          statically link the third-party .a archives; only system
#          frameworks (libSystem, libc++abi, ...) remain dynamic.
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
        # Linux portable static build (glibc-based).  Detailed flags are set
        # in the target_link_libraries section below.
        message(STATUS "STATIC_BUILD is ON (Linux/glibc) - producing a portable basevar (libstdc++/libgcc/htslib/compression libs all bundled, glibc dynamic)")
    endif()
endif()

# Ensure htslib is built before the main executable.
# For static builds on Linux:
#   * --disable-libcurl  : we only read local files, no remote URLs.
#   * --without-libdeflate: htslib's ./configure auto-detects libdeflate
#                          when it is present on the build host (e.g.
#                          GitHub-hosted ubuntu-22.04 runners have it
#                          preinstalled), and then emits libhts.a with
#                          libdeflate_* references.  Since we do not
#                          link -ldeflate (and want to keep the binary
#                          dependency-light), we explicitly disable it.
#                          htslib falls back to zlib transparently with
#                          no functional loss.
if(STATIC_BUILD AND NOT APPLE)
    set(HTSLIB_CONFIGURE_FLAGS --disable-libcurl --without-libdeflate)
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
    # ---- Linux portable partial-static build (glibc-based) ----
    # Build host: Ubuntu (CI uses ubuntu-22.04, glibc 2.35).
    # Strategy:
    #   * libstdc++ + libgcc are bundled via -static-libstdc++ -static-libgcc
    #   * htslib + zlib + bzip2 + lzma + openssl are linked as .a archives
    #   * glibc remains DYNAMIC (the binary still depends on the host's
    #     /lib64/ld-linux-x86-64.so.2 + libc.so.6, but those exist on every
    #     glibc-based Linux distribution).
    # This produces a binary that runs unmodified on any Linux with
    # glibc >= 2.31 (Ubuntu 20.04+, Debian 11+, RHEL/CentOS 8+, Fedora 32+).
    find_library(ZLIB_STATIC   NAMES libz.a     REQUIRED)
    find_library(BZ2_STATIC    NAMES libbz2.a   REQUIRED)
    find_library(LZMA_STATIC   NAMES liblzma.a  REQUIRED)
    find_library(SSL_STATIC    NAMES libssl.a   REQUIRED)
    find_library(CRYPTO_STATIC NAMES libcrypto.a REQUIRED)

    message(STATUS "Linux static libs:")
    message(STATUS "  zlib=${ZLIB_STATIC} bz2=${BZ2_STATIC} lzma=${LZMA_STATIC}")
    message(STATUS "  ssl=${SSL_STATIC} crypto=${CRYPTO_STATIC}")

    # Bundle the C++ runtime statically; glibc stays dynamic.
    set(CMAKE_EXE_LINKER_FLAGS
        "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc")

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
