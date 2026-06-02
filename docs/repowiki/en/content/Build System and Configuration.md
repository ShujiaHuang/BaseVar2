# Build System and Configuration

<cite>
**Referenced Files in This Document**
- [CMakeLists.txt](file://CMakeLists.txt)
- [README.md](file://README.md)
- [.github/workflows/build.yml](file://.github/workflows/build.yml)
- [.gitmodules](file://.gitmodules)
- [htslib/config.mk](file://htslib/config.mk)
- [htslib/htslib.mk](file://htslib/htslib.mk)
- [htslib/htslib_static.mk](file://htslib/htslib_static.mk)
- [htslib/configure.ac](file://htslib/configure.ac)
- [src/version.h.in](file://src/version.h.in)
- [src/version.h](file://src/version.h)
</cite>

## Update Summary
**Changes Made**
- Updated to reflect replacement of embedded htslib with external git submodule integration pointing to samtools/htslib repository
- Removed references to old embedded build system and replaced with proper submodule usage
- Updated installation instructions to reflect git submodule configuration
- Enhanced CI/CD workflow to handle recursive submodule checkout
- Updated platform compatibility requirements to glibc ≥ 2.35 for Ubuntu 22.04
- Simplified build system by leveraging external htslib submodule

## Table of Contents
1. [Introduction](#introduction)
2. [Project Structure](#project-structure)
3. [Core Components](#core-components)
4. [Architecture Overview](#architecture-overview)
5. [Detailed Component Analysis](#detailed-component-analysis)
6. [Static Build Configuration](#static-build-configuration)
7. [Dependency Analysis](#dependency-analysis)
8. [Performance Considerations](#performance-considerations)
9. [Troubleshooting Guide](#troubleshooting-guide)
10. [Conclusion](#conclusion)

## Introduction
This document explains BaseVar2's build system and configuration management with a focus on:
- CMake build configuration and control flow
- Dependency management for external htslib submodule and system libraries
- Platform-specific considerations (Linux vs macOS)
- Integration with external htslib submodule from samtools/htslib repository
- Compilation and linking requirements
- Runtime dependencies and development environment setup
- Compiler requirements and optimization flags
- **NEW**: Replaced embedded htslib with external git submodule integration for better maintainability and updates
- **NEW**: Enhanced CI/CD workflow with proper submodule handling for reliable builds
- **NEW**: Updated platform compatibility requirements to glibc ≥ 2.35 for Ubuntu 22.04

## Project Structure
BaseVar2 integrates a C++ application with an external htslib submodule from the samtools/htslib repository. The build system orchestrates:
- Version header generation via CMake configure_file
- Custom CMake target to build htslib using autotools from the submodule
- Inclusion of htslib headers and linking against the static htslib archive
- Platform-specific compiler flags and system libraries
- **NEW**: External submodule integration with proper git submodule configuration
- **NEW**: Enhanced CI/CD workflow with recursive submodule checkout support

```mermaid
graph TB
A["CMakeLists.txt"] --> B["Configure version header<br/>src/version.h.in → src/version.h<br/>(v2.2.3)"]
A --> C["Custom target Buildhts<br/>autoreconf/configure/make in htslib submodule"]
A --> D["Collect sources<br/>src/*.cpp + src/io/*.cpp"]
A --> E["Executable basevar"]
A --> F["Static Build Option<br/>STATIC_BUILD=ON"]
A --> G["Link libraries<br/>htslib static + system libs"]
A --> H["External htslib Submodule<br/>samtools/htslib repository"]
A --> I["Enhanced CI/CD Workflow<br/>Recursive submodule checkout"]
C --> J["htslib/libhts.a (from submodule)"]
F --> K["Platform-specific static flags"]
G --> J
G --> L["System libs<br/>pthread, z, bz2, m, lzma, curl, ssl, crypto"]
K --> M["Linux: -static-libstdc++ -static-libgcc + bundled libs"]
K --> N["macOS: static .a paths"]
H --> O["External repository integration<br/>Better maintainability"]
I --> P["Parallel dynamic and static builds<br/>with proper submodule handling"]
```

**Diagram sources**
- [CMakeLists.txt:17-20](file://CMakeLists.txt#L17-L20)
- [CMakeLists.txt:32-36](file://CMakeLists.txt#L32-L36)
- [CMakeLists.txt:52-55](file://CMakeLists.txt#L52-L55)
- [CMakeLists.txt:59](file://CMakeLists.txt#L59)
- [CMakeLists.txt:49](file://CMakeLists.txt#L49)
- [CMakeLists.txt:23](file://CMakeLists.txt#L23)
- [CMakeLists.txt:43-46](file://CMakeLists.txt#L43-L46)
- [CMakeLists.txt:158-192](file://CMakeLists.txt#L158-L192)
- [CMakeLists.txt:159-167](file://CMakeLists.txt#L159-L167)
- [.github/workflows/build.yml:32-35](file://.github/workflows/build.yml#L32-L35)
- [.github/workflows/build.yml:110-113](file://.github/workflows/build.yml#L110-L113)
- [.gitmodules:1-4](file://.gitmodules#L1-L4)

**Section sources**
- [CMakeLists.txt:1-197](file://CMakeLists.txt#L1-L197)
- [.github/workflows/build.yml:1-183](file://.github/workflows/build.yml#L1-L183)
- [.gitmodules:1-4](file://.gitmodules#L1-L4)

## Core Components
- CMake configuration enforces C++17, sets optimization flags, and configures platform-specific flags.
- A custom target builds htslib using autotools from the external submodule and places the resulting static library into the expected path.
- The main executable links against the static htslib archive and system libraries.
- Version metadata is templated into a generated header during configuration.
- **NEW**: External htslib submodule integration with proper git submodule configuration pointing to samtools/htslib repository.
- **NEW**: Enhanced CI/CD workflow with recursive submodule checkout support for reliable builds.
- **NEW**: Platform compatibility requirement updated to glibc ≥ 2.35 for Ubuntu 22.04.

Key implementation references:
- C++ standard and flags: [CMakeLists.txt:5](file://CMakeLists.txt#L5), [CMakeLists.txt:26-29](file://CMakeLists.txt#L26-L29)
- htslib build target: [CMakeLists.txt:32-36](file://CMakeLists.txt#L32-L36)
- Include directories and library linkage: [CMakeLists.txt:39-46](file://CMakeLists.txt#L39-L46), [CMakeLists.txt:49](file://CMakeLists.txt#L49), [CMakeLists.txt:61](file://CMakeLists.txt#L61)
- Version header generation: [CMakeLists.txt:17-20](file://CMakeLists.txt#L17-L20), [src/version.h.in:1-13](file://src/version.h.in#L1-L13), [src/version.h:1-13](file://src/version.h#L1-L13)
- Static build option: [CMakeLists.txt:23](file://CMakeLists.txt#L23), [CMakeLists.txt:58-75](file://CMakeLists.txt#L58-L75)
- **NEW**: External submodule configuration: [.gitmodules:1-4](file://.gitmodules#L1-L4)
- **NEW**: Enhanced CI/CD workflow: [.github/workflows/build.yml:32-35](file://.github/workflows/build.yml#L32-L35), [.github/workflows/build.yml:110-113](file://.github/workflows/build.yml#L110-L113)

**Section sources**
- [CMakeLists.txt:1-197](file://CMakeLists.txt#L1-L197)
- [src/version.h.in:1-13](file://src/version.h.in#L1-L13)
- [src/version.h:1-13](file://src/version.h#L1-L13)
- [.github/workflows/build.yml:32-35](file://.github/workflows/build.yml#L32-L35)
- [.github/workflows/build.yml:110-113](file://.github/workflows/build.yml#L110-L113)
- [.gitmodules:1-4](file://.gitmodules#L1-L4)

## Architecture Overview
The build architecture ties together CMake, external htslib submodule integration, and system libraries. The external submodule approach provides better maintainability and access to the latest htslib features while ensuring consistent builds across different environments. **NEW**: Enhanced CI/CD workflow with proper submodule handling ensures reliable builds with external dependencies.

```mermaid
sequenceDiagram
participant Dev as "Developer"
participant CMake as "CMakeLists.txt"
participant Git as "Git Submodule"
participant HTS as "External htslib (samtools/htslib)"
participant GCC as "Compiler/linker"
participant SYS as "System libs"
participant CI as "GitHub Actions"
participant Runner as "ubuntu-22.04 Runner"
participant Release as "GitHub Releases"
Dev->>CMake : Configure (C++17, flags, STATIC_BUILD?)
CMake->>Git : Checkout external htslib submodule
Git-->>CMake : htslib submodule (samtools/htslib)
CMake->>HTS : add_custom_target(Buildhts)<br/>autoreconf -i && ./configure && make
HTS-->>CMake : libhts.a (static)
CMake->>GCC : add_executable(basevar) with src/*.cpp + src/io/*.cpp
CMake->>GCC : add_dependencies(basevar Buildhts)
alt STATIC_BUILD=OFF
GCC->>HTS : Link libhts.a
GCC->>SYS : Link pthread, z, bz2, m, lzma, curl, [ssl, crypto]
else STATIC_BUILD=ON
GCC->>HTS : Link libhts.a with static flags
GCC->>SYS : Link system libs with bundled C++ runtime
note over GCC : Bundles libstdc++ + libgcc<br/>keeps glibc dynamic
end
GCC-->>Dev : bin/basevar
Note over CI : Enhanced CI/CD with recursive submodule<br/>checkout and glibc 2.35 compatibility
CI->>Runner : Build static binary (Linux)
Runner->>Runner : Build with -static-libstdc++ -static-libgcc
CI->>CI : Job 1 : Build & Test (dynamic linking)
CI->>CI : Job 2 : Static Build (ubuntu-22.04)
CI->>Release : Upload release assets with tag targeting
```

**Diagram sources**
- [CMakeLists.txt:32-36](file://CMakeLists.txt#L32-L36)
- [CMakeLists.txt:52-55](file://CMakeLists.txt#L52-L55)
- [CMakeLists.txt:59](file://CMakeLists.txt#L59)
- [CMakeLists.txt:61](file://CMakeLists.txt#L61)
- [CMakeLists.txt:43-46](file://CMakeLists.txt#L43-L46)
- [CMakeLists.txt:158-192](file://CMakeLists.txt#L158-L192)
- [.github/workflows/build.yml:32-35](file://.github/workflows/build.yml#L32-L35)
- [.github/workflows/build.yml:110-113](file://.github/workflows/build.yml#L110-L113)
- [.github/workflows/build.yml:82-93](file://.github/workflows/build.yml#L82-L93)
- [.github/workflows/build.yml:114-141](file://.github/workflows/build.yml#L114-141)
- [.github/workflows/build.yml:19-28](file://.github/workflows/build.yml#L19-L28)

## Detailed Component Analysis

### CMake Build Configuration
- Enforces C++17 and sets common compile flags including optimization level 3 and position-independent code.
- Adds a platform-specific linker flag for macOS to disable compact unwind information.
- Generates a version header from a template using CMake's configure_file mechanism.
- Declares a custom target to build htslib via autotools from the external submodule and ensures it completes before building the main executable.
- Includes htslib headers and links the static htslib archive with system libraries.
- **NEW**: External htslib submodule integration with proper git submodule configuration.
- **NEW**: Enhanced CI/CD workflow with recursive submodule checkout support.

```mermaid
flowchart TD
Start(["CMake configure"]) --> Std["Set C++17 standard"]
Std --> Flags["Set CXX flags:<br/>-O3 -fPIC (+ macOS -Wl,-no_compact_unwind)"]
Flags --> StaticCheck{"STATIC_BUILD?"}
StaticCheck --> |ON| StaticCfg["Apply platform-specific static flags"]
StaticCheck --> |OFF| VerHdr["configure_file(version.h.in → version.h)<br/>v2.2.3"]
StaticCfg --> LinuxStatic["Linux: -static-libstdc++ -static-libgcc<br/>+ bundled compression libs"]
StaticCfg --> MacStatic["macOS: static .a paths strategy"]
StaticCfg --> Ubuntu2204["ubuntu-22.04: glibc 2.35<br/>+ dynamic glibc support"]
StaticCfg --> VerHdr
Ubuntu2204 --> CurlFlags["Direct system curl linking"]
VerHdr --> HTSBuild["add_custom_target(Buildhts)<br/>autoreconf -i && ./configure && make<br/>(from external htslib submodule)"]
HTSBuild --> Srcs["Collect sources:<br/>src/*.cpp + src/io/*.cpp"]
Srcs --> Exec["add_executable(basevar)"]
Exec --> Deps["add_dependencies(basevar Buildhts)"]
Deps --> Link["target_link_libraries(basevar)<br/>libhts.a + system libs"]
Link --> End(["bin/basevar ready"])
```

**Diagram sources**
- [CMakeLists.txt:5](file://CMakeLists.txt#L5)
- [CMakeLists.txt:26-29](file://CMakeLists.txt#L26-L29)
- [CMakeLists.txt:17-20](file://CMakeLists.txt#L17-L20)
- [CMakeLists.txt:23](file://CMakeLists.txt#L23)
- [CMakeLists.txt:58-75](file://CMakeLists.txt#L58-L75)
- [CMakeLists.txt:158-192](file://CMakeLists.txt#L158-L192)
- [CMakeLists.txt:32-36](file://CMakeLists.txt#L32-L36)
- [CMakeLists.txt:52-55](file://CMakeLists.txt#L52-L55)
- [CMakeLists.txt:59](file://CMakeLists.txt#L59)
- [CMakeLists.txt:61](file://CMakeLists.txt#L61)

**Section sources**
- [CMakeLists.txt:1-197](file://CMakeLists.txt#L1-L197)

### External htslib Submodule Integration
- htslib is integrated as a git submodule pointing to the samtools/htslib repository, providing better maintainability and access to the latest features.
- The custom target builds htslib using autotools from the submodule and produces a static library at a fixed path expected by the main build.
- The main build includes htslib headers and links against the static archive.
- htslib's Makefile-based build supports optional plugins and compression libraries; the static archive is sufficient for BaseVar2's needs.
- **NEW**: External submodule integration with proper git submodule configuration pointing to samtools/htslib repository.
- **NEW**: Enhanced CI/CD workflow with recursive submodule checkout support for reliable builds.
- **NEW**: Improved maintainability and access to htslib updates without embedding code.

```mermaid
graph LR
A["CMakeLists.txt: Buildhts"] --> B["External htslib submodule<br/>(samtools/htslib)"]
B --> C["htslib/configure + make"]
C --> D["htslib/libhts.a"]
E["CMakeLists.txt: include_directories(htslib)"] --> F["Headers visible to basevar"]
G["CMakeLists.txt: target_link_libraries(... libhts.a ...)"] --> D
H["Enhanced CI/CD Workflow"] --> I["Recursive submodule checkout<br/>submodules: recursive"]
I --> G
```

**Diagram sources**
- [CMakeLists.txt:32-36](file://CMakeLists.txt#L32-L36)
- [CMakeLists.txt:39](file://CMakeLists.txt#L39)
- [CMakeLists.txt:49](file://CMakeLists.txt#L49)
- [CMakeLists.txt:61](file://CMakeLists.txt#L61)
- [.github/workflows/build.yml:32-35](file://.github/workflows/build.yml#L32-L35)
- [.github/workflows/build.yml:110-113](file://.github/workflows/build.yml#L110-L113)
- [.gitmodules:1-4](file://.gitmodules#L1-L4)

**Section sources**
- [CMakeLists.txt:32-36](file://CMakeLists.txt#L32-L36)
- [CMakeLists.txt:39](file://CMakeLists.txt#L39)
- [CMakeLists.txt:49](file://CMakeLists.txt#L49)
- [CMakeLists.txt:61](file://CMakeLists.txt#L61)
- [.github/workflows/build.yml:32-35](file://.github/workflows/build.yml#L32-L35)
- [.github/workflows/build.yml:110-113](file://.github/workflows/build.yml#L110-L113)
- [.gitmodules:1-4](file://.gitmodules#L1-L4)

### System Libraries and Platform-Specific Considerations
- Common system libraries: pthread, z (zlib), bz2 (bzip2), m (math), lzma (xz), curl.
- On non-macOS platforms, OpenSSL libraries (ssl and crypto) are appended.
- macOS-specific linker flag to disable compact unwind information is conditionally added.
- **NEW**: External htslib submodule integration provides consistent library availability across different environments.
- **NEW**: Linux static builds bundle libstdc++ and libgcc while keeping glibc dynamic for compatibility with glibc ≥ 2.35.
- **NEW**: Platform compatibility requirement updated to glibc ≥ 2.35 for modern Linux distribution support.
- **NEW**: Enhanced CI/CD workflow with proper submodule handling ensures reliable builds.

```mermaid
flowchart TD
A["Platform detection"] --> |macOS| M["Add macOS flag:<br/>-Wl,-no_compact_unwind"]
A --> |non-macOS| N["Append ssl + crypto"]
M --> L["Link system libs:<br/>pthread, z, bz2, m, lzma, curl, [ssl, crypto]"]
N --> L
A --> |STATIC_BUILD=ON| S["Static linking mode"]
S --> |Linux| LS["-static-libstdc++ -static-libgcc<br/>+ bundled compression libs"]
S --> |Linux ubuntu-22.04| LU["glibc 2.35 compatibility<br/>+ dynamic glibc support"]
S --> |macOS| MS["Static .a paths strategy"]
LS --> CurlDeps["Direct system curl linking"]
LU --> CurlDeps
LS --> SL["Portable static executable<br/>+ eliminates compatibility issues"]
LU --> LC["Ubuntu 22.04 base<br/>+ modern distribution support"]
MS --> SM["Best-effort static linking"]
```

**Diagram sources**
- [CMakeLists.txt:27-29](file://CMakeLists.txt#L27-L29)
- [CMakeLists.txt:44-46](file://CMakeLists.txt#L44-L46)
- [CMakeLists.txt:46-62](file://CMakeLists.txt#L46-L62)
- [CMakeLists.txt:158-192](file://CMakeLists.txt#L158-L192)
- [CMakeLists.txt:159-167](file://CMakeLists.txt#L159-L167)

**Section sources**
- [CMakeLists.txt:27-29](file://CMakeLists.txt#L27-L29)
- [CMakeLists.txt:44-46](file://CMakeLists.txt#L44-L46)
- [CMakeLists.txt:46-62](file://CMakeLists.txt#L46-L62)
- [CMakeLists.txt:158-192](file://CMakeLists.txt#L158-L192)
- [CMakeLists.txt:159-167](file://CMakeLists.txt#L159-L167)

### Version Header Generation
- A CMake configure_file step transforms a version template into a generated header containing project metadata.
- The template fields are populated from CMake project version and metadata variables.
- **NEW**: Version updated to 2.2.3 across all build configuration files.

```mermaid
sequenceDiagram
participant CMake as "CMake configure"
participant Tmpl as "src/version.h.in"
participant Out as "src/version.h"
CMake->>Tmpl : Read template
CMake->>Out : Write generated header with version info<br/>(v2.2.3)
```

**Diagram sources**
- [CMakeLists.txt:17-20](file://CMakeLists.txt#L17-L20)
- [src/version.h.in:1-13](file://src/version.h.in#L1-L13)
- [src/version.h:1-13](file://src/version.h#L1-L13)

**Section sources**
- [CMakeLists.txt:17-20](file://CMakeLists.txt#L17-L20)
- [src/version.h.in:1-13](file://src/version.h.in#L1-L13)
- [src/version.h:1-13](file://src/version.h#L1-L13)

### Enhanced CI/CD Build Workflow
- GitHub Actions builds on ubuntu-latest and macos-latest runners with enhanced support for external htslib submodule.
- **NEW**: Enhanced CI/CD workflow with proper submodule handling using `submodules: recursive` configuration.
- **NEW**: External htslib submodule integration ensures consistent library availability across different build environments.
- **NEW**: Platform compatibility requirement updated to glibc ≥ 2.35 for modern Linux distribution support.
- Installs dependencies via package managers (apt on Linux, Homebrew on macOS).
- Configures CMake with Release build type and enables testing.
- Builds and runs tests, then packages the binary artifact.
- **NEW**: Separate static build jobs for Linux (ubuntu-22.04) and macOS with improved compatibility.
- **NEW**: Enhanced verification process showing glibc version and symbol dependencies.

```mermaid
sequenceDiagram
participant GH as "GitHub Actions"
participant Perm as "Permissions Config"
participant Ubuntu as "ubuntu-latest runner"
participant Mac as "macos-latest runner"
participant Runner as "ubuntu-22.04 Runner"
participant Cfg as "CMake configure"
participant Bld as "CMake build"
participant Tst as "ctest"
GH->>Perm : Configure permissions : contents : write
GH->>Ubuntu : Job 1 : Build & Test (dynamic linking)
GH->>Mac : Job 1 : Build & Test (dynamic linking)
GH->>Runner : Job 2 : Static Build (Linux/ubuntu-22.04)
Runner->>Runner : Build with -static-libstdc++ -static-libgcc
Runner->>Runner : Verify glibc 2.35 compatibility
GH->>Cfg : cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
GH->>Bld : cmake --build build --config Release
GH->>Tst : ctest -C Release --output-on-failure
GH->>GH : Upload artifacts (bin/basevar)
GH->>GH : Upload Release Assets with tag targeting
```

**Diagram sources**
- [.github/workflows/build.yml:12-13](file://.github/workflows/build.yml#L12-L13)
- [.github/workflows/build.yml:19-28](file://.github/workflows/build.yml#L19-L28)
- [.github/workflows/build.yml:82-93](file://.github/workflows/build.yml#L82-L93)
- [.github/workflows/build.yml:114-141](file://.github/workflows/build.yml#L114-L141)
- [.github/workflows/build.yml:25-48](file://.github/workflows/build.yml#L25-L48)
- [.github/workflows/build.yml:50-52](file://.github/workflows/build.yml#L50-L52)
- [.github/workflows/build.yml:54-63](file://.github/workflows/build.yml#L54-L63)
- [.github/workflows/build.yml:177-183](file://.github/workflows/build.yml#L177-L183)

**Section sources**
- [.github/workflows/build.yml:1-183](file://.github/workflows/build.yml#L1-L183)

## Static Build Configuration

### STATIC_BUILD CMake Option
BaseVar2 introduces a comprehensive static linking capability through the `STATIC_BUILD` CMake option. This feature enables the creation of portable, self-contained binaries that can run on systems without the required shared libraries.

**Key Features:**
- **Cross-platform support**: Different strategies for Linux and macOS
- **NEW**: External htslib submodule integration provides consistent static library availability
- **NEW**: Enhanced CI/CD workflow with proper submodule handling for reliable builds
- **Best-effort static linking**: macOS static linking with system framework exceptions
- **NEW**: Platform compatibility requirement updated to glibc ≥ 2.35

### Linux Static Linking Strategy with Ubuntu 22.04
On Linux, the static build uses native Ubuntu 22.04 runner providing glibc 2.35 compatibility. The strategy bundles libstdc++ and libgcc while keeping glibc dynamic for maximum compatibility. **NEW**: External htslib submodule integration ensures consistent library availability across different host systems.

```cmake
if(STATIC_BUILD AND NOT APPLE)
    # ---- Linux portable partial-static build (glibc-based) ----
    # Build host: Ubuntu 22.04 (glibc 2.35).
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
endif()
```

**Benefits:**
- **NEW**: External htslib submodule integration ensures consistent library availability
- **NEW**: Compatible with Ubuntu 22.04+, Debian 12+, Fedora 36+, openSUSE Tumbleweed (glibc ≥ 2.35)
- **NEW**: Ensures compatibility with modern Linux distributions and enterprise environments
- **NEW**: Enhanced CI/CD workflow with proper submodule handling eliminates build inconsistencies
- **NEW**: Reduced runtime dependencies while maintaining broad compatibility
- **NEW**: Can run on minimal Docker containers with glibc ≥ 2.35

**Requirements:**
- **NEW**: Native Ubuntu 22.04 runner (glibc 2.35) for reliable static linking
- All required libraries must be available as static archives (zlib-static, bzip2-static, etc.)
- **NEW**: External htslib submodule integration ensures consistent library availability
- May limit certain features that rely on dynamic loading (DNS via NSS)
- **NEW**: glibc ≥ 2.35 requirement ensures compatibility with modern Linux distributions
- **NEW**: Enhanced CI/CD workflow eliminates Docker-based complexity

### macOS Static Linking Strategy
macOS has stricter limitations on static linking. The implementation uses a best-effort approach:

```cmake
if(STATIC_BUILD AND APPLE)
    # On macOS, Apple does NOT support fully static executables.
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
endif()
```

**Approach Details:**
- Static linking of compression libraries (zlib, bzip2, lzma) for version stability
- Dynamic linking of system frameworks (curl, ssl, crypto) due to Apple limitations
- Uses Homebrew paths for locating static archives
- Maintains compatibility with macOS system libraries

**Limitations:**
- Cannot create fully static executables on macOS due to Apple restrictions
- Some advanced features may require dynamic libraries
- Requires Homebrew installation of static library dependencies

### Static Build Verification
The enhanced GitHub Actions workflow includes improved verification steps for static binaries using native Ubuntu 22.04 runner. **NEW**: Enhanced verification process with glibc 2.35 compatibility and proper submodule handling.

**Linux Verification with Ubuntu 22.04:**
```bash
echo "=== Build environment ==="
gcc --version | head -1
g++ --version | head -1
cmake --version | head -1
ldd --version | head -1

echo "=== Verifying glibc 2.35 compatibility ==="
ldd --version | head -1

echo "=== Configuring (STATIC_BUILD=ON) ==="
cmake -B build-static -DSTATIC_BUILD=ON -DCMAKE_BUILD_TYPE=Release

echo "=== Building ==="
cmake --build build-static --config Release

echo "=== Verifying binary ==="
file bin/basevar
ldd bin/basevar || true
./bin/basevar | head -3
```

**macOS Verification:**
```bash
echo "=== Binary info ==="
file bin/basevar
echo "=== Dynamic library dependencies ==="
otool -L bin/basevar
```

**Section sources**
- [CMakeLists.txt:23](file://CMakeLists.txt#L23)
- [CMakeLists.txt:58-75](file://CMakeLists.txt#L58-L75)
- [CMakeLists.txt:126-156](file://CMakeLists.txt#L126-L156)
- [CMakeLists.txt:158-192](file://CMakeLists.txt#L158-L192)
- [.github/workflows/build.yml:82-93](file://.github/workflows/build.yml#L82-L93)
- [.github/workflows/build.yml:114-141](file://.github/workflows/build.yml#L114-L141)
- [.github/workflows/build.yml:158-162](file://.github/workflows/build.yml#L158-L162)

## Dependency Analysis
- Internal dependencies:
  - BaseVar executable depends on the htslib static archive from the external submodule.
  - The build ensures htslib is built before the executable.
- External dependencies:
  - System libraries: pthread, z, bz2, m, lzma, curl, ssl, crypto (conditional).
  - **NEW**: External htslib submodule integration with proper git submodule configuration.
  - **NEW**: Enhanced CI/CD workflow with recursive submodule checkout support.
  - **NEW**: Platform compatibility requirement updated to glibc ≥ 2.35.
- htslib internals:
  - htslib's Makefile exposes public headers and build targets for static/shared libraries and tools.
  - Optional plugins and compression libraries are selectable at configure time.
  - **NEW**: External htslib submodule integration ensures consistent static library generation.
  - **NEW**: Enhanced CI/CD workflow with proper submodule handling guarantees compatibility.

```mermaid
graph TB
subgraph "BaseVar2"
EXE["basevar (executable)"]
HDR["src headers"]
end
subgraph "External htslib Submodule"
LHTS["libhts.a (from samtools/htslib)"]
PUB["Public headers"]
SUBMODULE["git submodule<br/>samtools/htslib"]
END
subgraph "System libs"
P["pthread"]
Z["zlib"]
BZ["bzip2"]
M["m (libm)"]
LZ["lzma"]
CR["curl (system)"]
SSL["ssl"]
CRYPTO["crypto"]
end
EXE --> LHTS
EXE --> P
EXE --> Z
EXE --> BZ
EXE --> M
EXE --> LZ
EXE --> CR
EXE --> SSL
EXE --> CRYPTO
HDR --> PUB
SUBMODULE --> LHTS
subgraph "Enhanced CI/CD Workflow"
RECURSIVE["Recursive submodule checkout<br/>submodules: recursive"]
GLIBC["glibc 2.35 compatibility"]
COMPAT["Modern distribution support"]
end
RECURSIVE --> SUBMODULE
GLIBC --> COMPAT
COMPAT --> EXE
```

**Diagram sources**
- [CMakeLists.txt:49](file://CMakeLists.txt#L49)
- [CMakeLists.txt:61](file://CMakeLists.txt#L61)
- [CMakeLists.txt:43-46](file://CMakeLists.txt#L43-L46)
- [htslib/htslib.mk:57-83](file://htslib/htslib.mk#L57-L83)
- [CMakeLists.txt:158-192](file://CMakeLists.txt#L158-L192)
- [CMakeLists.txt:160-171](file://CMakeLists.txt#L160-L171)
- [.gitmodules:1-4](file://.gitmodules#L1-L4)
- [.github/workflows/build.yml:32-35](file://.github/workflows/build.yml#L32-L35)
- [.github/workflows/build.yml:110-113](file://.github/workflows/build.yml#L110-L113)

**Section sources**
- [CMakeLists.txt:43-61](file://CMakeLists.txt#L43-L61)
- [htslib/htslib.mk:57-83](file://htslib/htslib.mk#L57-L83)
- [htslib/htslib_static.mk:1-3](file://htslib/htslib_static.mk#L1-L3)
- [.gitmodules:1-4](file://.gitmodules#L1-L4)
- [.github/workflows/build.yml:32-35](file://.github/workflows/build.yml#L32-L35)
- [.github/workflows/build.yml:110-113](file://.github/workflows/build.yml#L110-L113)

## Performance Considerations
- Optimization flags:
  - The build uses aggressive optimization level 3 for performance-sensitive bioinformatics processing.
- Position-independent code:
  - -fPIC is enabled to support shared library compatibility and modern linkers.
- Platform-specific tuning:
  - macOS disables compact unwind information to avoid linker issues on some systems.
- Compression and networking:
  - htslib's autotools configuration detects CPU features and enables SIMD codecs when available, improving decompression performance.
- **NEW**: External htslib submodule integration reduces build complexity and improves maintainability.
- **NEW**: Native Ubuntu 22.04 static linking overhead:
  - Bundles only libstdc++ and libgcc while keeping glibc dynamic, reducing binary size.
  - **NEW**: Native Ubuntu 22.04 runner ensures consistent performance across different host systems.
  - May slightly reduce startup time by avoiding dynamic loading.
  - **NEW**: glibc 2.35 compatibility ensures compatibility with modern Linux distributions.
  - **NEW**: Enhanced CI/CD workflow eliminates Docker-based complexity and potential performance overhead.

Recommendations:
- Prefer Release builds for production use.
- Ensure sufficient virtual memory for multi-threaded variant calling; the project description indicates modest per-thread memory usage when tuned appropriately.
- **NEW**: Use static builds for environments with limited library availability or strict security requirements.
- **NEW**: For Linux deployments, consider native Ubuntu 22.04 runner for optimal static binary compatibility.
- **NEW**: The external htslib submodule approach eliminates build inconsistencies and improves maintainability.
- **NEW**: glibc 2.35 compatibility ensures compatibility with modern Linux distributions and enterprise environments.

**Section sources**
- [CMakeLists.txt:26](file://CMakeLists.txt#L26)
- [CMakeLists.txt:27-29](file://CMakeLists.txt#L27-L29)
- [htslib/configure.ac:86-191](file://htslib/configure.ac#L86-L191)

## Troubleshooting Guide

### Static Build Specific Issues

**Linux Static Build Failures with Ubuntu 22.04**
- **Symptom**: Linker errors indicating missing static archives or dependency issues
- **Resolution**: Ensure all required static libraries are available in Ubuntu 22.04 environment. The workflow includes proper installation of zlib1g-dev, libbz2-dev, liblzma-dev, libssl-dev packages.
- **References**: [.github/workflows/build.yml:118-121](file://.github/workflows/build.yml#L118-L121), [CMakeLists.txt:175-179](file://CMakeLists.txt#L175-L179)

**macOS Static Build Library Not Found**
- **Symptom**: CMake cannot locate static archives for zlib, bzip2, or lzma
- **Resolution**: Install dependencies via Homebrew and ensure they're built as static libraries. The build expects Homebrew paths. Note that curl, ssl, and crypto remain dynamic on macOS.
- **References**: [CMakeLists.txt:135-143](file://CMakeLists.txt#L135-L143)

**macOS Fully Static Limitation**
- **Symptom**: Static build still has dynamic dependencies
- **Resolution**: This is expected behavior on macOS. The build statically links third-party libraries while keeping system frameworks dynamic. Use the verification steps to confirm behavior.
- **References**: [CMakeLists.txt:59-69](file://CMakeLists.txt#L59-L69)

### Ubuntu 22.04 Compatibility Issues

**Symptom**: Static Linux binary incompatible with systems using older glibc versions
**Cause**: Platform compatibility requirement updated to glibc ≥ 2.35 (Ubuntu 22.04 runner)
**Solution**: Ensure deployment environment meets glibc ≥ 2.35 requirement (Ubuntu 22.04+, Debian 12+, Fedora 36+, openSUSE Tumbleweed)
**Verification**: Check glibc version with `ldd --version`

**References**: [CMakeLists.txt:159-167](file://CMakeLists.txt#L159-L167)

### CI/CD Workflow Issues

**GitHub Actions Permissions Error**
- **Symptom**: Release asset upload fails with permission denied
- **Resolution**: The workflow includes `permissions: contents: write` configuration. Ensure the workflow has the necessary permissions to upload release assets.
- **References**: [.github/workflows/build.yml:12-13](file://.github/workflows/build.yml#L12-L13)

**Static Build Job Dependencies**
- **Symptom**: Static build jobs wait for completion of dynamic build jobs
- **Resolution**: The static build jobs now run independently without `needs:` dependencies, enabling parallel execution. This improves CI/CD performance and reliability.
- **References**: [.github/workflows/build.yml:19-28](file://.github/workflows/build.yml#L19-L28)

**Release Asset Upload Issues**
- **Symptom**: Release assets not uploaded or incorrect tag targeting
- **Resolution**: The workflow uses explicit tag targeting with `${{ github.ref_name }}` for release assets. Verify that the release event is triggered with the correct tag name.
- **References**: [.github/workflows/build.yml:72-77](file://.github/workflows/build.yml#L72-L77), [.github/workflows/build.yml:177-183](file://.github/workflows/build.yml#L177-L183)

### External Submodule Integration Issues

**Submodule Checkout Failures**
- **Symptom**: htslib submodule not found or build fails
- **Resolution**: Ensure proper submodule initialization with `git submodule update --init --recursive`. The CI workflow uses `submodules: recursive` for automatic checkout.
- **References**: [.github/workflows/build.yml:32-35](file://.github/workflows/build.yml#L32-L35), [.github/workflows/build.yml:110-113](file://.github/workflows/build.yml#L110-L113)

**Submodule Update Issues**
- **Symptom**: Build fails after repository updates
- **Resolution**: Run `git submodule update --recursive` to sync with the latest htslib version. Development updates require submodule synchronization.
- **References**: [README.md:480-487](file://README.md#L480-L487)

**External Repository Access Issues**
- **Symptom**: Network connectivity problems with samtools/htslib repository
- **Resolution**: Ensure network access to GitHub. The submodule points to `git@github.com:samtools/htslib.git`. For HTTPS access, use `https://github.com/samtools/htslib.git`.
- **References**: [.gitmodules:1-4](file://.gitmodules#L1-L4)

### General Build Issues

**htslib autotools failures**
- Symptom: configure or make errors in htslib submodule.
- Resolution: The project's manual installation notes indicate that certain test-related errors can be ignored; the library still builds successfully. Retry autotools steps if transient network issues caused initial failure.
- References: [README.md:158](file://README.md#L158)

**Missing system libraries**
- **Symptom**: Linker errors for missing libraries (e.g., ssl, crypto, bz2, lzma, curl).
- **Resolution**: Install the required development packages using your OS package manager. The CI workflow demonstrates the exact packages installed on Linux and macOS, including Ubuntu 22.04 static packages.
- **References**: [.github/workflows/build.yml:36-41](file://.github/workflows/build.yml#L36-L41)

**macOS-specific linker issues**
- **Symptom**: Link-time errors related to compact unwind or symbol visibility.
- **Resolution**: The build adds a macOS-specific flag to disable compact unwind; ensure you are using the latest Xcode command line tools.
- **References**: [CMakeLists.txt:27-29](file://CMakeLists.txt#L27-L29)

**Manual linking differences**
- **Symptom**: Differences between automated CMake build and manual g++ invocation.
- **Resolution**: The manual script shows both approaches (using -lhts vs linking the static archive). Ensure include paths and library order match the CMake configuration.
- **References**: [README.md:142-156](file://README.md#L142-L156)

**Version header not updated**
- **Symptom**: Version metadata not reflecting project version.
- **Resolution**: Re-run CMake configure to regenerate the version header from the template. Version updated to 2.2.3 across all build configuration files.
- **References**: [CMakeLists.txt:17-20](file://CMakeLists.txt#L17-L20), [src/version.h.in:1-13](file://src/version.h.in#L1-L13), [src/version.h:1-13](file://src/version.h#L1-L13)

**Static Binary Verification Issues**
- **Symptom**: Static binary still shows dynamic dependencies
- **Resolution**: On macOS, this is expected behavior. The build statically links third-party libraries while keeping system frameworks dynamic. Use the verification steps in the CI workflow to confirm behavior.
- **References**: [.github/workflows/build.yml:158-162](file://.github/workflows/build.yml#L158-L162)

**Ubuntu 22.04 Container Issues**
- **Symptom**: Static build fails in Ubuntu 22.04 runner environment
- **Resolution**: Ensure all required static packages are installed (zlib1g-dev, libbz2-dev, liblzma-dev, libssl-dev). Clean stale htslib artifacts before building.
- **References**: [.github/workflows/build.yml:118-121](file://.github/workflows/build.yml#L118-L121)

**Ubuntu 22.04 Compatibility Issues**
- **Symptom**: Static Linux binary incompatible with older systems
- **Cause**: Platform compatibility requirement updated to glibc ≥ 2.35 (Ubuntu 22.04 runner)
- **Solution**: Deploy on systems meeting glibc ≥ 2.35 requirement (Ubuntu 22.04+, Debian 12+, Fedora 36+, openSUSE Tumbleweed)
- **Verification**: Check with `ldd --version` and ensure glibc version ≥ 2.35
- **References**: [CMakeLists.txt:159-167](file://CMakeLists.txt#L159-L167)

**Section sources**
- [.github/workflows/build.yml:12-13](file://.github/workflows/build.yml#L12-L13)
- [.github/workflows/build.yml:19-28](file://.github/workflows/build.yml#L19-L28)
- [.github/workflows/build.yml:72-77](file://.github/workflows/build.yml#L72-L77)
- [.github/workflows/build.yml:177-183](file://.github/workflows/build.yml#L177-L183)
- [CMakeLists.txt:135-143](file://CMakeLists.txt#L135-L143)
- [CMakeLists.txt:59-69](file://CMakeLists.txt#L59-L69)
- [README.md:158](file://README.md#L158)
- [.github/workflows/build.yml:36-41](file://.github/workflows/build.yml#L36-L41)
- [CMakeLists.txt:27-29](file://CMakeLists.txt#L27-L29)
- [README.md:142-156](file://README.md#L142-L156)
- [CMakeLists.txt:17-20](file://CMakeLists.txt#L17-L20)
- [src/version.h.in:1-13](file://src/version.h.in#L1-L13)
- [src/version.h:1-13](file://src/version.h#L1-L13)
- [.github/workflows/build.yml:158-162](file://.github/workflows/build.yml#L158-L162)
- [.github/workflows/build.yml:118-121](file://.github/workflows/build.yml#L118-L121)
- [.github/workflows/build.yml:138-141](file://.github/workflows/build.yml#L138-L141)
- [CMakeLists.txt:159-167](file://CMakeLists.txt#L159-L167)
- [.gitmodules:1-4](file://.gitmodules#L1-L4)
- [README.md:480-487](file://README.md#L480-L487)

## Conclusion
BaseVar2's build system centers on a streamlined CMake configuration that integrates external htslib via git submodule from samtools/htslib, enforces C++17, applies performance-oriented compiler flags, and links against essential system libraries. **NEW**: The introduction of external htslib submodule integration provides better maintainability, access to the latest htslib features, and consistent builds across different environments.

**NEW**: Enhanced CI/CD workflow with proper submodule handling ensures reliable builds with external dependencies. The new approach uses glibc 2.35 compatibility ensuring compatibility with modern Linux distributions while eliminating critical build inconsistencies that could arise from embedded htslib maintenance.

The enhanced CI/CD workflow represents a significant improvement in build automation and release management. **NEW**: The workflow now features recursive submodule checkout support, improved platform compatibility with glibc ≥ 2.35 requirement, and enhanced verification processes. These enhancements ensure consistent builds across different environments and eliminate the need for manual submodule management.

The external htslib submodule integration eliminates the need to maintain embedded code while providing access to the latest htslib features and bug fixes. The direct system curl linking and bundled C++ runtime libraries provide better compatibility while maintaining portability across modern Linux distributions with glibc ≥ 2.35.

The CI workflow validates builds on ubuntu-latest, macos-latest, and native Ubuntu 22.04 runners, ensuring broad platform compatibility and reliability. The enhanced GitHub Actions workflow now includes dedicated static build jobs with native Ubuntu 22.04 specifications and comprehensive verification steps, demonstrating the practical benefits of the external submodule approach for distribution and deployment scenarios.

Following the troubleshooting guidance and using the documented manual linking examples will resolve most build issues. For development, ensure the required system libraries are installed and use Release builds for optimal performance. **NEW**: For environments requiring zero-dependency deployments or strict security policies, utilize the STATIC_BUILD option with the appropriate platform-specific configuration, particularly leveraging native Ubuntu 22.04 runner for Linux static builds with its improved compatibility and elimination of build inconsistencies. The enhanced CI/CD workflow ensures reliable and efficient distribution of both dynamic and static binaries through parallel execution and improved release asset management.

**NEW**: Version 2.2.3 introduces critical stability improvements through the external htslib submodule integration, making BaseVar2 more robust and compatible across diverse deployment environments while eliminating the build inconsistencies present in the previous embedded approach. The glibc 2.35 compatibility ensures compatibility with modern Linux distributions and enterprise environments, while the enhanced CI/CD workflow provides improved reliability and performance through proper submodule handling.