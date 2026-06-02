# File Format I/O System

<cite>
**Referenced Files in This Document**
- [bam.h](file://src/io/bam.h)
- [bam.cpp](file://src/io/bam.cpp)
- [bam_header.h](file://src/io/bam_header.h)
- [bam_header.cpp](file://src/io/bam_header.cpp)
- [bam_record.h](file://src/io/bam_record.h)
- [bam_record.cpp](file://src/io/bam_record.cpp)
- [fasta.h](file://src/io/fasta.h)
- [fasta.cpp](file://src/io/fasta.cpp)
- [vcf.h](file://src/io/vcf.h)
- [vcf.cpp](file://src/io/vcf.cpp)
- [vcf_header.h](file://src/io/vcf_header.h)
- [vcf_header.cpp](file://src/io/vcf_header.cpp)
- [vcf_record.h](file://src/io/vcf_record.h)
- [vcf_record.cpp](file://src/io/vcf_record.cpp)
- [hts_utils.h](file://src/io/hts_utils.h)
- [utils.h](file://src/io/utils.h)
- [utils.cpp](file://src/io/utils.cpp)
- [CMakeLists.txt](file://CMakeLists.txt)
- [.gitmodules](file://.gitmodules)
- [README.md](file://README.md)
</cite>

## Update Summary
**Changes Made**
- Enhanced documentation for htslib submodule integration and dependency management
- Updated file format handling documentation to reflect new external htslib submodule approach
- Added comprehensive coverage of static vs dynamic linking strategies
- Expanded troubleshooting guidance for htslib integration issues
- Updated installation and build instructions reflecting submodule management

## Table of Contents
1. [Introduction](#introduction)
2. [Project Structure](#project-structure)
3. [Core Components](#core-components)
4. [Architecture Overview](#architecture-overview)
5. [Detailed Component Analysis](#detailed-component-analysis)
6. [htslib Submodule Integration](#htslib-submodule-integration)
7. [Dependency Analysis](#dependency-analysis)
8. [Performance Considerations](#performance-considerations)
9. [Troubleshooting Guide](#troubleshooting-guide)
10. [Conclusion](#conclusion)

## Introduction
This document describes BaseVar2's file format I/O system with emphasis on:
- BAM/CRAM alignment file processing
- FASTA reference genome handling
- VCF output generation

It explains how the system integrates with htslib for maximum compatibility and performance, details stream processing and compression support, and documents the internal data structures used for file representation and conversion between formats. Error handling, validation procedures, and performance optimization techniques for large file processing are covered.

**Updated** Enhanced documentation reflects improved integration of external htslib submodule with comprehensive dependency management approach, including both static and dynamic linking strategies for optimal deployment flexibility.

## Project Structure
The I/O system is organized into focused modules under src/io:
- BAM/CRAM: [bam.h](file://src/io/bam.h), [bam.cpp](file://src/io/bam.cpp), [bam_header.h](file://src/io/bam_header.h), [bam_header.cpp](file://src/io/bam_header.cpp), [bam_record.h](file://src/io/bam_record.h), [bam_record.cpp](file://src/io/bam_record.cpp)
- FASTA: [fasta.h](file://src/io/fasta.h), [fasta.cpp](file://src/io/fasta.cpp)
- VCF/BCF: [vcf.h](file://src/io/vcf.h), [vcf.cpp](file://src/io/vcf.cpp), [vcf_header.h](file://src/io/vcf_header.h), [vcf_header.cpp](file://src/io/vcf_header.cpp), [vcf_record.h](file://src/io/vcf_record.h), [vcf_record.cpp](file://src/io/vcf_record.cpp)
- Utilities: [hts_utils.h](file://src/io/hts_utils.h), [utils.h](file://src/io/utils.h), [utils.cpp](file://src/io/utils.cpp)

```mermaid
graph TB
subgraph "I/O Layer"
BAM["BAM/CRAM<br/>bam.h/cpp"]
FASTA["FASTA<br/>fasta.h/cpp"]
VCF["VCF/BCF<br/>vcf.h/cpp"]
end
subgraph "Data Models"
BH["BamHeader<br/>bam_header.h/cpp"]
BR["BamRecord<br/>bam_record.h/cpp"]
VH["VCFHeader<br/>vcf_header.h/cpp"]
VR["VCFRecord<br/>vcf_record.h/cpp"]
end
subgraph "Utilities"
HU["hts_utils.h"]
U["utils.h/cpp"]
end
BAM --> BH
BAM --> BR
VCF --> VH
VCF --> VR
BAM --> HU
VCF --> HU
FASTA --> HU
BAM --> U
VCF --> U
FASTA --> U
```

**Diagram sources**
- [bam.h:20-149](file://src/io/bam.h#L20-L149)
- [bam.cpp:1-167](file://src/io/bam.cpp#L1-L167)
- [bam_header.h:18-121](file://src/io/bam_header.h#L18-L121)
- [bam_header.cpp:1-102](file://src/io/bam_header.cpp#L1-L102)
- [bam_record.h:49-455](file://src/io/bam_record.h#L49-L455)
- [bam_record.cpp:1-551](file://src/io/bam_record.cpp#L1-L551)
- [fasta.h:14-96](file://src/io/fasta.h#L14-L96)
- [fasta.cpp:1-122](file://src/io/fasta.cpp#L1-L122)
- [vcf.h:29-184](file://src/io/vcf.h#L29-L184)
- [vcf.cpp:1-227](file://src/io/vcf.cpp#L1-L227)
- [vcf_header.h:31-242](file://src/io/vcf_header.h#L31-L242)
- [vcf_header.cpp:1-275](file://src/io/vcf_header.cpp#L1-L275)
- [vcf_record.h:31-525](file://src/io/vcf_record.h#L31-L525)
- [vcf_record.cpp:1-800](file://src/io/vcf_record.cpp#L1-L800)
- [hts_utils.h:17-61](file://src/io/hts_utils.h#L17-L61)
- [utils.h:19-205](file://src/io/utils.h#L19-L205)
- [utils.cpp:1-142](file://src/io/utils.cpp#L1-L142)

**Section sources**
- [bam.h:1-149](file://src/io/bam.h#L1-L149)
- [vcf.h:1-184](file://src/io/vcf.h#L1-L184)
- [fasta.h:1-96](file://src/io/fasta.h#L1-L96)

## Core Components
- BAM/CRAM reader/writer:
  - [Bam:23-145](file://src/io/bam.h#L23-L145): Opens, indexes, iterates, and streams SAM/BAM/CRAM records via htslib. Supports region queries and CRAM with reference.
  - [BamHeader:22-118](file://src/io/bam_header.h#L22-L118): Wraps sam_hdr_t, exposes contig names/lengths, sample name extraction.
  - [BamRecord:49-455](file://src/io/bam_record.h#L49-L455): Encapsulates bam1_t, provides flags, mapping info, CIGAR parsing, aligned pairs, tags, and utilities.

- FASTA reader:
  - [Fasta:16-91](file://src/io/fasta.h#L16-L91): Indexes FASTA via faidx, supports region fetch and basic metadata.

- VCF/BCF reader/writer:
  - [VCFFile:29-179](file://src/io/vcf.h#L29-L179): Opens, reads/writes VCF/BCF, manages header and iterators, supports region queries.
  - [VCFHeader:31-242](file://src/io/vcf_header.h#L31-L242): Manages bcf_hdr_t with shared ownership, sample/contig access, header manipulation.
  - [VCFRecord:31-525](file://src/io/vcf_record.h#L31-L525): Manages bcf1_t with shared ownership, accessors/mutators for core fields, INFO/FILTER/FORMAT.

- Utilities:
  - [hts_utils.h:17-61](file://src/io/hts_utils.h#L17-L61): Format detection helpers and CRAM detection.
  - [utils.h:21-205](file://src/io/utils.h#L21-L205), [utils.cpp:1-142](file://src/io/utils.cpp#L1-L142): File/path helpers, string conversions, splitting, and I/O utilities.

**Section sources**
- [bam.h:20-149](file://src/io/bam.h#L20-L149)
- [bam_header.h:18-121](file://src/io/bam_header.h#L18-L121)
- [bam_record.h:49-455](file://src/io/bam_record.h#L49-L455)
- [fasta.h:14-96](file://src/io/fasta.h#L14-L96)
- [vcf.h:29-184](file://src/io/vcf.h#L29-L184)
- [vcf_header.h:31-242](file://src/io/vcf_header.h#L31-L242)
- [vcf_record.h:31-525](file://src/io/vcf_record.h#L31-L525)
- [hts_utils.h:17-61](file://src/io/hts_utils.h#L17-L61)
- [utils.h:19-205](file://src/io/utils.h#L19-L205)
- [utils.cpp:1-142](file://src/io/utils.cpp#L1-L142)

## Architecture Overview
The system builds on htslib for native support of SAM/BAM/CRAM and VCF/BCF, with thin C++ wrappers that:
- Manage resource lifecycles (RAII) using smart pointers for header/record objects
- Provide region-based streaming via htslib iterators
- Support transparent compression (BGZF/GZIP) and indexing (.bai/.csi/.crai/.tbi)
- Offer convenience APIs for common operations (e.g., fetching FASTA sequences by region)

```mermaid
sequenceDiagram
participant App as "Application"
participant VCF as "VCFFile"
participant Hdr as "VCFHeader"
participant Rec as "VCFRecord"
participant Hts as "htslib"
App->>VCF : open("file.vcf.gz", "r")
VCF->>Hts : hts_open/read header
Hts-->>VCF : bcf_hdr_t*
VCF->>Hdr : wrap header
App->>VCF : index_load()
VCF->>Hts : hts_idx_load2()
Hts-->>VCF : hts_idx_t*
App->>VCF : fetch(region)
VCF->>Hts : bcf_itr_querys()
Hts-->>VCF : hts_itr_t*
loop iterate
App->>VCF : read(Rec)
VCF->>Hts : bcf_itr_next()/bcf_read()
Hts-->>VCF : status
VCF-->>App : Rec
end
```

**Diagram sources**
- [vcf.cpp:8-57](file://src/io/vcf.cpp#L8-L57)
- [vcf.cpp:87-107](file://src/io/vcf.cpp#L87-L107)
- [vcf.cpp:109-161](file://src/io/vcf.cpp#L109-L161)
- [vcf.cpp:164-198](file://src/io/vcf.cpp#L164-L198)

**Section sources**
- [vcf.h:29-179](file://src/io/vcf.h#L29-L179)
- [vcf.cpp:1-227](file://src/io/vcf.cpp#L1-L227)

## Detailed Component Analysis

### BAM/CRAM Processing
- File opening and mode handling:
  - Supports reading/writing modes compatible with htslib, including binary/text, compression, and CRAM with reference.
  - CRAM requires a reference path; the wrapper sets the FAI filename via htslib.
- Indexing and iteration:
  - Loads .bai/.csi as needed; creates iterators for region queries.
  - Iterators support string regions and coordinate-based queries.
- Streaming and record access:
  - Sequential or iterator-based reading via sam_read1/sam_itr_next.
  - Records expose flags, mapping info, CIGAR, query sequences/qualities, and tags.

```mermaid
classDiagram
class Bam {
-string _fname
-string _mode
-int _io_status
-string _reference_path
-samFile* _fp
-hts_idx_t* _idx
-hts_itr_t* _itr
-BamHeader _hdr
+Bam(fn, mode, ref_fn)
+index_build(min_shift)
+index_load()
+fetch(region) bool
+fetch(seq_id, beg, end) bool
+read(BamRecord) int
+io_status() int
}
class BamHeader {
-sam_hdr_t* _h
+BamHeader(samFile*)
+BamHeader(fn, ref_fn)
+name2id(name) int
+get_sample_name() string
+write(samFile*) int
}
class BamRecord {
-bam1_t* _b
+load_read(samFile*, sam_hdr_t*) int
+next_read(samFile*, hts_itr_t*) int
+query_sequence() string
+query_qual(offset) string
+cigar() string
+get_cigar_blocks() vector
+get_alignment_blocks() vector
+get_aligned_pairs(ref) vector
+has_tag(tag) bool
+get_tag(tag) string
}
Bam --> BamHeader : "owns header"
Bam --> BamRecord : "reads/writes"
```

**Diagram sources**
- [bam.h:23-145](file://src/io/bam.h#L23-L145)
- [bam_header.h:22-118](file://src/io/bam_header.h#L22-L118)
- [bam_record.h:49-455](file://src/io/bam_record.h#L49-L455)

**Section sources**
- [bam.h:20-149](file://src/io/bam.h#L20-L149)
- [bam.cpp:6-46](file://src/io/bam.cpp#L6-L46)
- [bam_header.h:18-118](file://src/io/bam_header.h#L18-L118)
- [bam_header.cpp:5-39](file://src/io/bam_header.cpp#L5-L39)
- [bam_record.h:49-455](file://src/io/bam_record.h#L49-L455)
- [bam_record.cpp:90-120](file://src/io/bam_record.cpp#L90-L120)

### FASTA Reference Handling
- Index-based access:
  - Uses faidx to load .fai indices; supports bgzip-compressed FASTA.
  - Provides region fetch by string ("chr:start-end") or by coordinates.
- Thread safety:
  - The implementation notes that direct fetch is not thread-safe; use per-thread instances or guard access.

```mermaid
flowchart TD
Start(["Open FASTA"]) --> CheckIndex["Load .fai index"]
CheckIndex --> |Success| Ready["Ready for fetch"]
CheckIndex --> |Failure| Error["Throw invalid_argument"]
Ready --> FetchStr["Fetch by region string"]
Ready --> FetchCoord["Fetch by chr/start/end"]
FetchStr --> Validate["Validate region"]
FetchCoord --> Validate
Validate --> |OK| ReturnSeq["Return subsequence"]
Validate --> |Invalid| ThrowErr["Throw invalid_argument"]
```

**Diagram sources**
- [fasta.cpp:9-22](file://src/io/fasta.cpp#L9-L22)
- [fasta.cpp:50-95](file://src/io/fasta.cpp#L50-L95)

**Section sources**
- [fasta.h:14-96](file://src/io/fasta.h#L14-L96)
- [fasta.cpp:9-122](file://src/io/fasta.cpp#L9-L122)

### VCF/BCF Output Generation
- Reader/writer lifecycle:
  - Opening in read mode reads the header; in write mode writes the header immediately.
  - Index loading supports .tbi/.csi; region queries supported via iterators.
- Record model:
  - VCFRecord wraps bcf1_t with shared ownership; provides accessors for CHROM/POS/ID/REF/ALT/QUAL/FILTER/INFO/FORMAT.
  - Supports unpacking for efficient field access and mutation APIs for INFO/FORMAT updates.
- Header model:
  - VCFHeader wraps bcf_hdr_t with shared ownership; supports sample/contig access, subset operations, and header line manipulation.

```mermaid
classDiagram
class VCFFile {
-string _fname
-string _mode
-int _io_status
-htsFile* _fp
-hts_idx_t* _idx
-hts_itr_t* _itr
-VCFHeader _hdr
+VCFFile(fn, mode)
+VCFFile(fn, hdr, mode)
+index_load()
+fetch(region) bool
+fetch(seq_id, beg, end) bool
+read(VCFRecord) int
+write(VCFRecord) int
+io_status() int
}
class VCFHeader {
-shared_ptr<bcf_hdr_t> _hdr
+n_samples() int
+sample_names() vector<string>
+subset_samples(names) VCFHeader
+seq_name(id) string
+seq_id(name) int
+to_string() string
}
class VCFRecord {
-shared_ptr<bcf1_t> _b
+rid(hdr) int32_t
+chrom(hdr) string
+pos() hts_pos_t
+ref() string
+alt() vector<string>
+qual() float
+unpack(which) int
+get_info_* / get_format_* / update_* ...
+subset_samples(hdr, names) VCFRecord
}
VCFFile --> VCFHeader : "owns header"
VCFFile --> VCFRecord : "reads/writes"
```

**Diagram sources**
- [vcf.h:29-179](file://src/io/vcf.h#L29-L179)
- [vcf.cpp:8-57](file://src/io/vcf.cpp#L8-L57)
- [vcf_header.h:31-242](file://src/io/vcf_header.h#L31-L242)
- [vcf_record.h:31-525](file://src/io/vcf_record.h#L31-L525)

**Section sources**
- [vcf.h:29-184](file://src/io/vcf.h#L29-L184)
- [vcf.cpp:1-227](file://src/io/vcf.cpp#L1-L227)
- [vcf_header.h:31-242](file://src/io/vcf_header.h#L31-L242)
- [vcf_header.cpp:1-275](file://src/io/vcf_header.cpp#L1-L275)
- [vcf_record.h:31-525](file://src/io/vcf_record.h#L31-L525)
- [vcf_record.cpp:1-800](file://src/io/vcf_record.cpp#L1-L800)

## htslib Submodule Integration

### External Submodule Management
BaseVar2 utilizes htslib as an external Git submodule for enhanced dependency management and build flexibility. The integration follows these key approaches:

- **Submodule Declaration**: The htslib repository is declared as a submodule in [.gitmodules:1-4](file://.gitmodules#L1-L4) with the URL pointing to the official samtools/htslib repository.
- **Build Automation**: CMake handles htslib compilation through a custom target that executes autoreconf, configure, and make steps within the htslib directory.
- **Flexible Linking**: The build system supports both static and dynamic linking strategies depending on deployment requirements.

### Static vs Dynamic Linking Strategies

#### Static Build (Recommended for Portability)
The static build strategy bundles all dependencies including htslib and compression libraries into a single executable:

```mermaid
flowchart TD
A["Static Build Request"] --> B["Configure htslib with --disable-libcurl --without-libdeflate"]
B --> C["Compile htslib to libhts.a"]
C --> D["Link libhts.a statically"]
D --> E["Bundle zlib, bzip2, lzma, openssl statically"]
E --> F["Produce portable executable"]
```

**Diagram sources**
- [CMakeLists.txt:77-99](file://CMakeLists.txt#L77-L99)
- [CMakeLists.txt:111-113](file://CMakeLists.txt#L111-L113)

Key characteristics of static builds:
- **Linux**: Partial-static approach bundling libstdc++, libgcc, htslib, and compression libs; glibc remains dynamic for broad compatibility.
- **macOS**: Best-effort static linking with system frameworks remaining dynamic due to platform limitations.
- **Portability**: Executables run on various Linux distributions without external dependencies.

#### Dynamic Build (Development Flexibility)
The dynamic build links against system-installed htslib libraries:

```mermaid
flowchart LR
A["Dynamic Build"] --> B["Link against system libhts.so"]
B --> C["Runtime dependency resolution"]
C --> D["Requires htslib installation"]
D --> E["Development and testing flexibility"]
```

### Build Configuration and Dependencies

#### System Requirements
- **C++17 Compiler**: GCC 7+ or Apple Clang 10+
- **CMake**: Version 3.12 or higher
- **System Libraries**: zlib, bzip2, xz-utils, libcurl
- **Optional**: OpenSSL for certain features

#### Installation Methods
The project supports multiple installation approaches:

1. **Pre-built Binaries**: Download static binaries from GitHub Releases for immediate use
2. **Source Compilation**: Clone with submodules and build using CMake
3. **Manual Compilation**: Build htslib separately, then compile BaseVar2

**Section sources**
- [.gitmodules:1-4](file://.gitmodules#L1-L4)
- [CMakeLists.txt:77-99](file://CMakeLists.txt#L77-L99)
- [CMakeLists.txt:111-113](file://CMakeLists.txt#L111-L113)
- [CMakeLists.txt:194-196](file://CMakeLists.txt#L194-L196)
- [README.md:96-142](file://README.md#L96-L142)

## Dependency Analysis
- Internal dependencies:
  - BAM depends on BamHeader and BamRecord; both rely on htslib structs.
  - VCFFile depends on VCFHeader and VCFRecord; both rely on htslib structs.
  - FASTA depends on faidx from htslib.
- External dependencies:
  - htslib for SAM/BAM/CRAM, VCF/BCF, FAI indexing, and compression.
  - Standard C++ filesystem and iostream facilities for path and I/O utilities.
- **Updated** Enhanced dependency management through external htslib submodule with flexible linking strategies.

```mermaid
graph LR
HU["hts_utils.h"] --> BAM["Bam"]
HU --> VCF["VCFFile"]
HU --> FASTA["Fasta"]
Utils["utils.h/cpp"] --> BAM
Utils --> VCF
Utils --> FASTA
BAM --> BH["BamHeader"]
BAM --> BR["BamRecord"]
VCF --> VH["VCFHeader"]
VCF --> VR["VCFRecord"]
subgraph "htslib Integration"
HTS["libhts.a/.so"] --> HU
HTS --> BAM
HTS --> VCF
HTS --> FASTA
end
```

**Diagram sources**
- [hts_utils.h:17-61](file://src/io/hts_utils.h#L17-L61)
- [utils.h:19-205](file://src/io/utils.h#L19-L205)
- [utils.cpp:1-142](file://src/io/utils.cpp#L1-L142)
- [bam.h:13-18](file://src/io/bam.h#L13-L18)
- [vcf.h:13-19](file://src/io/vcf.h#L13-L19)
- [fasta.h:11-11](file://src/io/fasta.h#L11-L11)

**Section sources**
- [hts_utils.h:17-61](file://src/io/hts_utils.h#L17-L61)
- [utils.h:19-205](file://src/io/utils.h#L19-L205)
- [utils.cpp:1-142](file://src/io/utils.cpp#L1-L142)
- [bam.h:13-18](file://src/io/bam.h#L13-L18)
- [vcf.h:13-19](file://src/io/vcf.h#L13-L19)
- [fasta.h:11-11](file://src/io/fasta.h#L11-L11)

## Performance Considerations
- Streaming and indexing:
  - Use region-based queries with pre-built indices (.bai/.csi/.crai/.tbi) to avoid full-file scans.
  - Prefer BGZF-compressed inputs for SAM/BAM/VCF/BCF to enable random access via indices.
- Memory efficiency:
  - Header and record objects use shared ownership to minimize duplication and simplify cleanup.
  - Iterator-based reading avoids loading entire files into memory.
- Compression:
  - Transparent BGZF/GZIP support via htslib; choose appropriate compression levels for output writers.
- Parallelization:
  - The BAM module comments indicate care for thread-safety with iterators; avoid sharing a single file's iterator across threads. Instead, process separate files concurrently.
- **Updated** Enhanced performance through optimized htslib integration with static linking for reduced overhead and improved portability.

## Troubleshooting Guide
- File opening failures:
  - Verify readability and existence of input files; exceptions are thrown for missing files or open failures.
  - For CRAM, ensure the reference FASTA path is provided and readable.
- Index-related errors:
  - Ensure index files (.bai/.csi/.crai/.tbi/.tbi) exist and match the input file; the reader throws if index loading fails.
- Region queries:
  - Confirm region strings are valid and contig names exist in the header; invalid regions cause query failures.
- FASTA fetch:
  - Ensure the FASTA file is indexed; fetching by region throws on invalid inputs or empty results.
- Header validation:
  - For VCF writers, ensure the header is valid before opening; invalid headers lead to write failures.
- **Updated** htslib integration issues:
  - Verify htslib submodule initialization: `git submodule update --init --recursive`
  - Check build configuration for static vs dynamic linking requirements
  - Ensure proper linking of compression libraries (zlib, bzip2, lzma, openssl)
  - Validate system library availability for dynamic builds

**Section sources**
- [bam.cpp:6-46](file://src/io/bam.cpp#L6-L46)
- [bam.cpp:86-97](file://src/io/bam.cpp#L86-L97)
- [bam.cpp:103-135](file://src/io/bam.cpp#L103-L135)
- [fasta.cpp:9-22](file://src/io/fasta.cpp#L9-L22)
- [fasta.cpp:50-95](file://src/io/fasta.cpp#L50-L95)
- [vcf.cpp:8-57](file://src/io/vcf.cpp#L8-L57)
- [vcf.cpp:87-107](file://src/io/vcf.cpp#L87-L107)
- [vcf.cpp:109-161](file://src/io/vcf.cpp#L109-L161)
- [CMakeLists.txt:77-99](file://CMakeLists.txt#L77-L99)
- [README.md:100-108](file://README.md#L100-L108)

## Conclusion
BaseVar2's I/O system leverages htslib to deliver robust, high-performance processing of alignment and variant data. The modular design with RAII-managed header and record objects, region-based streaming, and transparent compression support enables scalable workflows for large-scale genomics data. The FASTA interface simplifies reference access, while the VCF/BCF stack provides flexible reading and writing with strong validation and error reporting.

**Updated** The enhanced integration with external htslib submodule provides improved dependency management, flexible linking strategies, and better deployment options. The static linking approach ensures portability across diverse environments, while dynamic linking offers development flexibility. This dual-strategy approach addresses both production deployment needs and development workflow requirements, making BaseVar2 more adaptable to various use cases and environments.