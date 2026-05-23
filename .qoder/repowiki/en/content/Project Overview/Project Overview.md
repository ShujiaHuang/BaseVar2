# Project Overview

<cite>
**Referenced Files in This Document**
- [README.md](file://README.md)
- [src/main.cpp](file://src/main.cpp)
- [src/version.h](file://src/version.h)
- [src/variant_caller.h](file://src/variant_caller.h)
- [src/basetype.h](file://src/basetype.h)
- [src/algorithm.h](file://src/algorithm.h)
- [scripts/create_pipeline.py](file://scripts/create_pipeline.py)
- [update_note.md](file://update_note.md)
- [.github/workflows/build.yml](file://.github/workflows/build.yml)
- [tests/data/sample_group.info](file://tests/data/sample_group.info)
- [CMakeLists.txt](file://CMakeLists.txt)
- [docs/assets/images/basevar_logo.svg](file://docs/assets/images/basevar_logo.svg)
- [docs/assets/images/basevar_Image.png](file://docs/assets/images/basevar_Image.png)
- [docs/assets/images/old_basevar_logo.svg](file://docs/assets/images/old_basevar_logo.svg)
</cite>

## Update Summary
**Changes Made**
- Enhanced visual presentation with larger logo (300px height) and improved branding
- Updated documentation to reflect the new visual identity system with SVG logo scaling
- Added documentation for the comprehensive visual asset management system
- Updated project structure section to reflect enhanced visual presentation capabilities

## Table of Contents
1. [Introduction](#introduction)
2. [Project Structure](#project-structure)
3. [Core Components](#core-components)
4. [Architecture Overview](#architecture-overview)
5. [Detailed Component Analysis](#detailed-component-analysis)
6. [Dependency Analysis](#dependency-analysis)
7. [Performance Considerations](#performance-considerations)
8. [Troubleshooting Guide](#troubleshooting-guide)
9. [Conclusion](#conclusion)
10. [Appendices](#appendices)

## Introduction
BaseVar2 is a specialized ultra-low-depth whole genome sequencing (ULDS/WGS) variant caller designed for non-invasive prenatal testing (NIPT) and related human genetic research. The latest 2.2.3 release maintains the high-performance C++ implementation while introducing a simplified build approach and honest system requirements. It targets single-nucleotide polymorphism (SNP) and insertion–deletion (Indel) detection from sub-single-read coverage data (<1x), enabling cost-effective, population-scale NIPT studies. The tool emphasizes:

- Accurate variant detection via likelihood-based inference
- Population-level allele frequency estimation from sparse coverage
- High-performance C++ implementation with substantial speed and memory improvements over the original Python version
- **Enhanced**: Improved visual branding with scalable SVG logo system featuring larger 300px height presentation for optimal landing page experience

**Updated** The project now features an enhanced visual identity system with a high-quality SVG logo that provides superior scalability and rendering quality compared to traditional raster formats. The logo is prominently displayed at 300px height in the README.md landing page, ensuring crisp presentation on high-DPI displays, mobile devices, and various screen sizes while maintaining excellent print quality. This enhanced visual presentation improves the project's professional appearance and user experience on the main landing page.

Key scientific and technical goals:
- Robust detection of rare variants in extremely shallow coverage
- Reliable estimation of minor allele frequencies (MAF) and genotype likelihoods
- Scalable, parallelized processing suitable for large cohorts and whole-genome analyses
- **Enhanced**: Honest system requirements with glibc >= 2.35 for Linux static binaries
- **Enhanced**: Modernized visual presentation with scalable branding assets optimized for documentation and web interfaces

Citation and publication details:
- Liu et al. (2024). Utilizing non-invasive prenatal test sequencing data for human genetic investigation. Cell Genomics 4(10), 100669. https://doi.org/10.1016/j.xgen.2024.100669

## Project Structure
High-level organization:
- CLI entrypoint and command routing
- Core variant caller engine and batching logic
- Mathematical/statistical modules for likelihoods, tests, and EM-based allele frequency updates
- I/O wrappers for FASTA, BAM/CRAM, and VCF
- Pipeline generation script for chromosome-wise or region-wise parallelization
- HTSlib integration for efficient NGS file parsing
- Build automation and CI workflows with simplified cross-platform support
- **Enhanced**: Scalable SVG logo system with optimized 300px height presentation for landing pages

```mermaid
graph TB
A["CLI Entrypoint<br/>src/main.cpp"] --> B["Caller Engine<br/>src/variant_caller.h"]
B --> C["Core Model: BaseType<br/>src/basetype.h"]
B --> D["Statistical Algorithms<br/>src/algorithm.h"]
B --> E["I/O Adapters<br/>src/io/*"]
B --> F["Population Grouping<br/>tests/data/sample_group.info"]
A --> G["Concat/Subset Tools<br/>src/concat.h, src/vcf_subset_samples.h"]
A --> H["Pipeline Generator<br/>scripts/create_pipeline.py"]
B --> I["HTSlib Integration<br/>htslib/*"]
J[".github/workflows/build.yml"] --> K["Enhanced Build & Test Automation"]
L["CMakeLists.txt<br/>Version 2.2.3"] --> M["Simplified Build System"]
N["SVG Logo System<br/>docs/assets/images/basevar_logo.svg"] --> O["300px Height Landing Page Presentation"]
P["PNG Fallback<br/>docs/assets/images/basevar_Image.png"] --> Q["Raster Compatibility"]
R["Legacy SVG<br/>docs/assets/images/old_basevar_logo.svg"] --> S["Backward Compatibility"]
```

**Diagram sources**
- [src/main.cpp:1-105](file://src/main.cpp#L1-L105)
- [src/variant_caller.h:1-180](file://src/variant_caller.h#L1-L180)
- [src/basetype.h:1-146](file://src/basetype.h#L1-L146)
- [src/algorithm.h:1-180](file://src/algorithm.h#L1-L180)
- [scripts/create_pipeline.py:1-103](file://scripts/create_pipeline.py#L1-L103)
- [.github/workflows/build.yml:1-183](file://.github/workflows/build.yml#L1-L183)
- [CMakeLists.txt:1-197](file://CMakeLists.txt#L1-L197)
- [docs/assets/images/basevar_logo.svg:1-100](file://docs/assets/images/basevar_logo.svg#L1-L100)
- [docs/assets/images/basevar_Image.png:1-100](file://docs/assets/images/basevar_Image.png#L1-L100)
- [docs/assets/images/old_basevar_logo.svg:1-100](file://docs/assets/images/old_basevar_logo.svg#L1-L100)

**Section sources**
- [README.md:1-488](file://README.md#L1-L488)
- [src/main.cpp:1-105](file://src/main.cpp#L1-L105)
- [scripts/create_pipeline.py:1-103](file://scripts/create_pipeline.py#L1-L103)
- [.github/workflows/build.yml:1-183](file://.github/workflows/build.yml#L1-L183)
- [CMakeLists.txt:1-197](file://CMakeLists.txt#L1-L197)
- [docs/assets/images/basevar_logo.svg:1-100](file://docs/assets/images/basevar_logo.svg#L1-L100)
- [docs/assets/images/basevar_Image.png:1-100](file://docs/assets/images/basevar_Image.png#L1-L100)
- [docs/assets/images/old_basevar_logo.svg:1-100](file://docs/assets/images/old_basevar_logo.svg#L1-L100)

## Core Components
- CLI and Commands
  - Subcommands: caller, concat, subsam
  - Version and author metadata exposed at runtime
- Variant Caller Engine
  - Region-aware batching and parallelization
  - Sample grouping support for population-specific allele frequency calculations
  - Output VCF generation with standardized headers and fields
- Statistical Foundation
  - Likelihood ratio testing (LRT) thresholds and quality filters
  - Fisher's exact test and chi-squared-based tests for significance
  - Expectation–Maximization (EM) algorithm for allele frequency estimation
- I/O Layer
  - FASTA reference handling
  - BAM/CRAM alignment parsing and filtering by MAPQ/Quality
  - BGZF-compressed VCF writing and header management
- Pipeline Utilities
  - Automated chromosome-wise or region-wise job splitting
  - Optional population group file for stratified analyses
- **Enhanced**: Scalable Visual Identity System
  - SVG logo with crisp rendering at any resolution, prominently displayed at 300px height
  - PNG fallback for compatibility with legacy systems and documentation builds
  - Legacy SVG support for backward compatibility
  - Consistent branding across documentation, presentations, and web interfaces with optimized landing page presentation

Practical examples (usage patterns):
- Single-sample or multi-sample ULDS calling with region selection and population grouping
- Batched processing via a sample list file
- Whole-genome parallelization using the pipeline generator
- **Enhanced**: Professional visual presentation in research publications and presentations with optimized logo scaling

**Section sources**
- [src/main.cpp:18-32](file://src/main.cpp#L18-L32)
- [src/version.h:1-13](file://src/version.h#L1-L13)
- [src/variant_caller.h:41-174](file://src/variant_caller.h#L41-L174)
- [src/algorithm.h:90-178](file://src/algorithm.h#L90-L178)
- [scripts/create_pipeline.py:26-94](file://scripts/create_pipeline.py#L26-L94)
- [tests/data/sample_group.info:1-44](file://tests/data/sample_group.info#L1-L44)
- [docs/assets/images/basevar_logo.svg:1-100](file://docs/assets/images/basevar_logo.svg#L1-L100)

## Architecture Overview
The system orchestrates a pipeline from input alignments to population-aware variant calls:
- Parse CLI and initialize caller parameters
- Load reference FASTA and optionally population groups
- Partition genome into regions and create batches
- For each batch:
  - Fetch aligned bases per sample and filter by quality/mapping criteria
  - Compute per-site likelihoods and perform LRT-based variant scoring
  - Estimate population-level allele frequencies via EM
  - Write VCF records with FORMAT and INFO fields
- Optionally concatenate outputs and subset samples post-run
- **Enhanced**: Consistent visual branding throughout the workflow presentation with optimized logo display

```mermaid
sequenceDiagram
participant U as "User"
participant CLI as "CLI (main.cpp)"
participant ENG as "Caller Engine (variant_caller.h)"
participant MOD as "Model (basetype.h)"
participant ALG as "Algorithms (algorithm.h)"
participant IO as "I/O (FASTA/BAM/VCF)"
participant PIPE as "Pipeline (create_pipeline.py)"
participant BRAND as "Visual Branding (300px Logo)"
U->>CLI : "basevar caller ..."
CLI->>ENG : "Parse args, init"
ENG->>IO : "Load reference, samples, regions"
ENG->>PIPE : "Split regions (optional)"
loop For each batch
ENG->>IO : "Fetch bases (MAPQ/Q filtering)"
ENG->>MOD : "Compute likelihoods"
MOD->>ALG : "LRT, EM updates"
ALG-->>MOD : "Posterior, MAF"
ENG->>IO : "Write VCF record"
end
CLI-->>U : "Report completion and timings"
U->>BRAND : "Display enhanced branded output"
```

**Diagram sources**
- [src/main.cpp:45-104](file://src/main.cpp#L45-L104)
- [src/variant_caller.h:120-137](file://src/variant_caller.h#L120-L137)
- [src/basetype.h:95-143](file://src/basetype.h#L95-L143)
- [src/algorithm.h:150-178](file://src/algorithm.h#L150-L178)
- [scripts/create_pipeline.py:75-94](file://scripts/create_pipeline.py#L75-L94)
- [docs/assets/images/basevar_logo.svg:1-100](file://docs/assets/images/basevar_logo.svg#L1-L100)

## Detailed Component Analysis

### CLI and Command Routing
- Provides top-level command dispatch (caller, concat, subsam)
- Prints version and author metadata
- Captures start/end timestamps and CPU time for performance reporting

```mermaid
flowchart TD
Start(["Process Start"]) --> CheckArgs{"Args present?"}
CheckArgs --> |No| Help["Print usage and exit"]
CheckArgs --> |Yes| ParseCmd["Parse command: caller/concat/subsam"]
ParseCmd --> Caller{"caller?"}
ParseCmd --> Concat{"concat?"}
ParseCmd --> Subsam{"subsam?"}
Caller --> RunCaller["Run variant caller"]
Concat --> RunConcat["Run concatenation"]
Subsam --> RunSubsam["Run subsample extraction"]
RunCaller --> End(["Exit"])
RunConcat --> End
RunSubsam --> End
```

**Diagram sources**
- [src/main.cpp:45-104](file://src/main.cpp#L45-L104)

**Section sources**
- [src/main.cpp:18-32](file://src/main.cpp#L18-L32)
- [src/main.cpp:45-104](file://src/main.cpp#L45-L104)
- [src/version.h:1-13](file://src/version.h#L1-L13)

### Variant Calling Engine
- Loads calling intervals and sample IDs
- Supports population grouping for stratified MAF computation
- Implements batch creation and per-region processing
- Integrates quality filters (min MAPQ, min base quality) and region selection

```mermaid
classDiagram
class BaseTypeRunner {
+run() int
-_make_calling_interval()
-_get_sample_id_from_bam()
-_get_popgroup_info()
-_create_batchfiles(...)
-_fetch_base_in_region(...)
-_seek_position(...)
-_write_record_to_batchfile(...)
-_variants_discovery(...)
-_variant_calling_unit(...)
-_basevar_caller(...)
-_basetype_caller_unit(...)
}
```

**Diagram sources**
- [src/variant_caller.h:41-174](file://src/variant_caller.h#L41-L174)

**Section sources**
- [src/variant_caller.h:41-174](file://src/variant_caller.h#L41-L174)

### Mathematical Foundation and Model
- LRT-based variant scoring with predefined thresholds
- Quality-based conversion constants for Phred-scale to log-space
- Fisher's exact test and chi-squared tests for significance
- EM algorithm for iterative allele frequency updates across samples

```mermaid
flowchart TD
A["Observed Bases per Site"] --> B["Compute Allele Likelihoods"]
B --> C["Likelihood Ratio Test (LRT)"]
C --> D{"Exceeds Threshold?"}
D --> |Yes| E["Call Variant + Qual Score"]
D --> |No| F["Skip/Mark Low-Quality"]
E --> G["Estimate Allele Frequencies (EM)"]
G --> H["Update Posterior Probabilities"]
H --> I["Write VCF Record"]
```

**Diagram sources**
- [src/basetype.h:29-143](file://src/basetype.h#L29-L143)
- [src/algorithm.h:150-178](file://src/algorithm.h#L150-L178)

**Section sources**
- [src/basetype.h:25-28](file://src/basetype.h#L25-L28)
- [src/algorithm.h:100-138](file://src/algorithm.h#L100-L138)
- [src/algorithm.h:150-178](file://src/algorithm.h#L150-L178)

### Pipeline Generation and Parallelization
- Splits a reference FASTA index into fixed-size genomic windows
- Generates shell commands to run caller on each window in parallel
- Supports population grouping and region filtering

```mermaid
sequenceDiagram
participant User as "User"
participant Script as "create_pipeline.py"
participant FS as "Filesystem"
participant Exec as "basevar caller"
User->>Script : "Provide ref.fai, delta, chrom, threads"
Script->>FS : "Read ref.fai"
loop For each chromosome and window
Script->>Exec : "Generate caller command with -r region"
Exec-->>FS : "Output VCF per window"
end
Script-->>User : "Shell script to execute"
```

**Diagram sources**
- [scripts/create_pipeline.py:26-94](file://scripts/create_pipeline.py#L26-L94)

**Section sources**
- [scripts/create_pipeline.py:26-94](file://scripts/create_pipeline.py#L26-L94)

### Enhanced Visual Identity System
**Updated** The project now features a comprehensive visual branding system designed for modern bioinformatics applications with enhanced presentation:

- **Primary SVG Logo**: High-resolution vector logo located at `docs/assets/images/basevar_logo.svg` featuring prominent 300px height display for optimal landing page presentation and crisp rendering at any scale
- **PNG Fallback**: Compatible raster backup at `docs/assets/images/basevar_Image.png` ensuring broad compatibility across documentation platforms and legacy systems
- **Legacy Support**: Backward compatibility maintained through `docs/assets/images/old_basevar_logo.svg` for systems requiring older logo formats
- **Scalability Benefits**: Vector graphics maintain quality across different display densities, from standard monitors to high-DPI screens and mobile devices
- **Cross-Platform Consistency**: Ensures uniform appearance across documentation websites, presentations, academic posters, and research materials
- **Optimized Landing Page Experience**: The 300px height logo provides professional prominence while maintaining fast loading times and excellent visual impact
- **Mobile Optimization**: Scales appropriately for mobile device displays and touch interfaces without quality degradation
- **Print Quality**: Maintains sharp detail in printed materials and academic publications regardless of scaling

The visual identity system supports:
- Research presentations with professional appearance and optimal visual impact
- Documentation websites with responsive design and enhanced user experience
- Academic posters and conference materials with consistent branding
- Social media and promotional content with scalable vector graphics
- Mobile applications and web interfaces with adaptive scaling capabilities

**Section sources**
- [docs/assets/images/basevar_logo.svg:1-100](file://docs/assets/images/basevar_logo.svg#L1-L100)
- [docs/assets/images/basevar_Image.png:1-100](file://docs/assets/images/basevar_Image.png#L1-L100)
- [docs/assets/images/old_basevar_logo.svg:1-100](file://docs/assets/images/old_basevar_logo.svg#L1-L100)

### Practical Examples and Target Use Cases
- NIPT cohort-wide SNP detection from ULDS data
- Population-level MAF estimation by grouping samples (e.g., regional or ethnic groups)
- Whole-genome analysis via automated chromosome-wise pipelines
- Downstream compatibility with tools expecting standard VCF FORMAT/INFO fields
- **Enhanced**: Professional visual presentation in research communications with optimized logo display

Target scenarios:
- Large-scale maternal plasma ULDS studies requiring speed and memory efficiency
- Cost-effective screening protocols leveraging shallow coverage
- Research applications needing robust variant detection and frequency estimation
- **Enhanced**: High-quality visual presentation for academic and commercial contexts with professional logo branding

**Section sources**
- [README.md:13-488](file://README.md#L13-L488)
- [tests/data/sample_group.info:1-44](file://tests/data/sample_group.info#L1-L44)

## Dependency Analysis
Internal dependencies:
- CLI depends on caller engine and I/O modules
- Caller engine depends on model (BaseType), algorithms, and I/O adapters
- Algorithms module provides shared statistical primitives
- Pipeline script depends on CLI binary availability
- **Enhanced**: Visual assets depend on proper asset management and deployment with optimized file sizes

External dependencies:
- HTSlib for NGS file parsing and compression
- Standard system libraries (pthread, zlib, bz2, lzma, curl, OpenSSL on non-macOS)
- **Enhanced**: Asset delivery systems for scalable visual content with optimized caching

```mermaid
graph LR
CLI["CLI (main.cpp)"] --> CALLER["Caller (variant_caller.h)"]
CALLER --> MODEL["Model (basetype.h)"]
CALLER --> ALG["Algorithms (algorithm.h)"]
CALLER --> IO["I/O (io/*)"]
CALLER --> PIPE["Pipeline (create_pipeline.py)"]
CALLER --> HTS["HTSlib"]
ALG --> HTS
ASSETS["Visual Assets<br/>SVG/PNG/OLD_LOGOS"] --> BRANDING["Enhanced Brand Presentation (300px)"]
BRANDING --> DOCS["Documentation"]
BRANDING --> WEBSITE["Web Interface"]
```

**Diagram sources**
- [src/main.cpp:12-16](file://src/main.cpp#L12-L16)
- [src/variant_caller.h:23-31](file://src/variant_caller.h#L23-L31)
- [src/algorithm.h:22](file://src/algorithm.h#L22)
- [docs/assets/images/basevar_logo.svg:1-100](file://docs/assets/images/basevar_logo.svg#L1-L100)

**Section sources**
- [CMakeLists.txt:31-62](file://CMakeLists.txt#L31-L62)
- [.github/workflows/build.yml:29-40](file://.github/workflows/build.yml#L29-L40)

## Performance Considerations
- C++ implementation delivers >10x speedup versus the Python predecessor with significantly reduced memory footprint
- Per-thread memory consumption scales with batch size and region size; tuning batch count and threads improves throughput
- Quality filters (min MAPQ and base quality) reduce noise and improve accuracy without heavy computational overhead
- Parallelization via region-wise pipelines accelerates whole-genome analyses
- **Updated**: Simplified build system reduces compilation time and build failures across platforms
- **Enhanced**: Optimized visual asset delivery for faster documentation loading and better user experience with 300px logo scaling

[No sources needed since this section provides general guidance]

## Troubleshooting Guide
Common issues and resolutions:
- Compilation failures due to missing system libraries (zlib, bz2, lzma, curl, OpenSSL on Linux)
  - Ensure all build dependencies are installed before configuring and building
- HTSlib submodule configuration errors
  - Re-run autotools and configure inside the htslib directory; warnings in tests are often benign
- Memory usage spikes
  - Reduce batch size or thread count; adjust batch count and region granularity
- Incorrect population grouping
  - Verify the population group file format matches sample IDs in BAM headers
- **Updated**: Static binary compatibility issues
  - Use pre-built static binaries for zero-dependency deployment
  - Linux static binaries require glibc >= 2.35 (Ubuntu 22.04+, Debian 12+, Fedora 36+)
  - macOS static binaries require macOS 12+ and have minimal system dependencies
- **Enhanced**: Visual asset loading issues
  - Ensure SVG and PNG assets are properly deployed in documentation builds with correct file paths
  - Verify asset dimensions match the expected 300px height specification for landing page display
  - Check browser compatibility for SVG rendering in older environments
  - Confirm fallback PNG assets are accessible when SVG rendering fails

**Section sources**
- [.github/workflows/build.yml:29-40](file://.github/workflows/build.yml#L29-L40)
- [README.md:57-84](file://README.md#L57-L84)
- [update_note.md:7-31](file://update_note.md#L7-L31)

## Conclusion
BaseVar2 advances the state-of-the-art for ULDS variant discovery by combining rigorous statistical modeling with high-performance C++ implementation. The 2.2.3 release maintains significant improvements in build reliability and cross-platform compatibility while adopting a simplified build approach. Its focus on accurate variant detection and reliable allele frequency estimation from <1x data makes it especially suited for NIPT and population-scale studies. The modular architecture, robust I/O layer, and automated parallelization pipeline enable scalable, reproducible workflows from shallow coverage datasets.

**Enhanced** The project's new SVG-based visual identity system with 300px height logo presentation provides professional, scalable branding that enhances the overall presentation quality for research communications, academic publications, and commercial applications. The combination of technical excellence and enhanced visual polish positions BaseVar2 as a comprehensive solution for modern bioinformatics workflows with optimal landing page experience and consistent brand identity across all platforms.

[No sources needed since this section summarizes without analyzing specific files]

## Appendices

### Citation and Publication Details
- Liu et al. (2024). Utilizing non-invasive prenatal test sequencing data for human genetic investigation. Cell Genomics 4(10), 100669. https://doi.org/10.1016/j.xgen.2024.100669

**Section sources**
- [README.md:13-26](file://README.md#L13-L26)

### Ultra-Low-Depth Sequencing Challenges (Beginner-Friendly)
- Coverage limitations: Very low read depth increases uncertainty in genotyping and reduces power to detect rare variants
- Base quality and mapping quality filtering: Essential to mitigate sequencing errors and alignment artifacts
- Allele dropout and bias: Population-level statistics require careful modeling to avoid false positives
- Computational efficiency: Scalable algorithms and parallelization are critical for whole-genome analyses

[No sources needed since this section provides general guidance]

### Version 2.2.3 Build Improvements
**Simplified Build System Features:**
- **Streamlined Static Linking**: Simplified approach focusing on glibc >= 2.35 compatibility
- **Honest System Requirements**: Clear glibc version requirements for Linux static binaries
- **Reduced Build Complexity**: Eliminated manylinux2014 compatibility claims in favor of realistic requirements
- **Improved Reliability**: Better error handling and clearer build instructions
- **Enhanced Visual Assets**: New SVG logo system with 300px height optimization for landing page presentation

**Static Binary Availability:**
- **Linux (x86_64)**: Partial-static binary with glibc >= 2.35 requirement
- **macOS (arm64/Intel)**: Best-effort static binary with minimal system requirements

**Enhanced Visual Identity Features:**
- **SVG Logo**: High-resolution vector logo with 300px height for optimal landing page display
- **PNG Fallback**: Compatible raster backup for broad compatibility across documentation platforms
- **Legacy Support**: Backward compatibility maintained for older systems and applications
- **Cross-Platform Consistency**: Uniform appearance across different devices, displays, and print media
- **Mobile Optimization**: Scales appropriately for mobile and tablet interfaces with crisp rendering
- **Professional Presentation**: Enhanced visual impact for research communications and academic contexts

**Section sources**
- [README.md:19-43](file://README.md#L19-L43)
- [.github/workflows/build.yml:82-183](file://.github/workflows/build.yml#L82-L183)
- [CMakeLists.txt:46-63](file://CMakeLists.txt#L46-L63)
- [docs/assets/images/basevar_logo.svg:1-100](file://docs/assets/images/basevar_logo.svg#L1-L100)
- [docs/assets/images/basevar_Image.png:1-100](file://docs/assets/images/basevar_Image.png#L1-L100)
- [docs/assets/images/old_basevar_logo.svg:1-100](file://docs/assets/images/old_basevar_logo.svg#L1-L100)

### Visual Asset Management
**Enhanced Asset Structure:**
- Primary SVG logo: `docs/assets/images/basevar_logo.svg` with optimized 300px height for landing page presentation
- Fallback PNG image: `docs/assets/images/basevar_Image.png` for compatibility with legacy systems and documentation builds
- Legacy SVG support: `docs/assets/images/old_basevar_logo.svg` for backward compatibility
- Asset dimensions optimized for different use cases and display densities
- Cross-platform compatibility ensured across modern browsers and documentation systems

**Usage Guidelines:**
- SVG preferred for scalable applications and high-DPI displays with optimal landing page presentation
- PNG fallback for legacy systems and environments without SVG support or documentation platform compatibility
- Consistent sizing recommendations for different presentation contexts including 300px height for main landing page
- Proper attribution and licensing considerations for public use across research and commercial applications

**Technical Specifications:**
- SVG vector format with embedded fonts and scalable elements optimized for web delivery
- PNG raster format with appropriate resolution for web and print with backward compatibility
- File size optimization for fast loading and distribution across documentation platforms
- Cross-browser compatibility verified across major platforms including modern documentation systems
- Responsive design principles applied for optimal display across different screen sizes and devices

**Section sources**
- [docs/assets/images/basevar_logo.svg:1-100](file://docs/assets/images/basevar_logo.svg#L1-L100)
- [docs/assets/images/basevar_Image.png:1-100](file://docs/assets/images/basevar_Image.png#L1-L100)
- [docs/assets/images/old_basevar_logo.svg:1-100](file://docs/assets/images/old_basevar_logo.svg#L1-L100)