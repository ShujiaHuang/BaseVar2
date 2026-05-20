# Citation and Publication Information

<cite>
**Referenced Files in This Document**
- [README.md](file://README.md)
- [CMakeLists.txt](file://CMakeLists.txt)
- [src/version.h](file://src/version.h)
- [src/main.cpp](file://src/main.cpp)
- [src/variant_caller.h](file://src/variant_caller.h)
- [src/algorithm.h](file://src/algorithm.h)
- [src/basetype.cpp](file://src/basetype.cpp)
</cite>

## Table of Contents
1. [Introduction](#introduction)
2. [Primary Publication Reference](#primary-publication-reference)
3. [Research Significance of NIPT Data](#research-significance-of-nipt-data)
4. [Citing BaseVar2 in Academic Works](#citing-basevar2-in-academic-works)
5. [Collaborative Research Team and Affiliations](#collaborative-research-team-and-affiliations)
6. [Methodology Acknowledgment Guidelines](#methodology-acknowledgment-guidelines)
7. [Technical Implementation Notes](#technical-implementation-notes)
8. [Conclusion](#conclusion)

## Introduction

This document provides comprehensive citation and publication information for BaseVar2 research, focusing on the primary publication in Cell Genomics (2024) by Huang et al. and associated guidelines for proper academic attribution. BaseVar2 is a specialized tool designed for variant calling using ultra-low-depth (<1x) sequencing data, particularly targeting non-invasive prenatal test (NIPT) data in human genetic studies.

## Primary Publication Reference

### Complete Citation Details

**Primary Reference:**
- Liu, S., Liu, Y., Gu, Y., Lin, X., Zhu, H., Liu, H., Xu, Z., Cheng, S., Lan, X., Li, L., Huang, M., Li, H., Nielsen, R., Davies, RW., Albrechtsen, A., Chen, GB., Qiu, X., Jin, X., **Huang, S.**, (2024). Utilizing non-invasive prenatal test sequencing data for human genetic investigation. *Cell Genomics* 4(10), 100669 [doi:10.1016/j.xgen.2024.100669](https://www.cell.com/cell-genomics/fulltext/S2666-979X(24)00288-X)

### Journal Information

- **Journal:** Cell Genomics
- **Volume:** 4
- **Issue:** 10
- **Article Number:** 100669
- **Year:** 2024
- **DOI:** 10.1016/j.xgen.2024.100669
- **Publisher:** Elsevier

### Mathematical Foundation Reference

The research builds upon foundational mathematical work referenced in the project documentation:
- [Mathematical explanations for maximum likelihood and likelihood ratio models](https://doi.org/10.1016/j.cell.2018.08.016)

**Section sources**
- [README.md:17](file://README.md#L17)

## Research Significance of NIPT Data

### Non-Invasive Prenatal Testing (NIPT) Applications

BaseVar2's utilization of NIPT sequencing data represents a significant advancement in human genetic investigation methodology. NIPT data offers several key advantages:

- **Non-invasive sampling:** Eliminates risks associated with invasive procedures like amniocentesis
- **Early detection capabilities:** Enables genetic screening during early pregnancy stages
- **Cost-effective screening:** Provides broad genomic coverage at reduced costs compared to traditional methods
- **Population-scale studies:** Facilitates large-scale genetic epidemiology research

### Scientific Impact

The research demonstrates how ultra-low-depth sequencing data can be effectively utilized for:
- Accurate polymorphism identification at genomic positions
- Reliable allele frequency calculations
- Population genetic studies using maternal plasma-derived DNA
- Development of computational tools for clinical applications

**Section sources**
- [README.md:9](file://README.md#L9)

## Citing BaseVar2 in Academic Works

### Primary Citation Requirement

When using BaseVar2 in published works, research papers, or academic presentations, cite the following primary reference:

**Huang, S., et al. (2024). Utilizing non-invasive prenatal test sequencing data for human genetic investigation. *Cell Genomics* 4(10), 100669.**

### Version-Specific Citations

For specific versions of BaseVar2, include the version number in your citation:

**Example citation format:**
- BaseVar2 v2.0.0 (2024). Software tool for variant calling from ultra-low-depth WGS data. Available at https://github.com/ShujiaHuang/basevar2

### Command-Line Tool Citation

When reporting results obtained through BaseVar2, cite both the tool and methodology:

**Recommended citation structure:**
- "Variants were called using BaseVar2 v2.0.0 (Huang et al., 2024) with default parameters for ultra-low-depth WGS data analysis."

### Academic Presentation Guidelines

For conference presentations and posters:
- Include the citation in presentation slides
- Provide software availability information
- Reference the specific version used
- Acknowledge funding sources if applicable

**Section sources**
- [README.md:15-17](file://README.md#L15-L17)
- [src/version.h:4-10](file://src/version.h#L4-L10)
- [CMakeLists.txt:9](file://CMakeLists.txt#L9)

## Collaborative Research Team and Affiliations

### Core Research Team Members

The BaseVar2 research involved an international collaborative team with diverse institutional affiliations:

**Primary Authors:**
- Liu, S., Liu, Y., Gu, Y., Lin, X., Zhu, H., Liu, H., Xu, Z., Cheng, S., Lan, X., Li, L., Huang, M., Li, H., Nielsen, R., Davies, RW., Albrechtsen, A., Chen, GB., Qiu, X., Jin, X., **Huang, S.**

### Institutional Affiliations

Based on the publication, the research team includes affiliations from multiple institutions across different countries, representing a truly international collaboration in computational genomics and bioinformatics research.

### Collaborative Research Impact

The collaborative nature of this research demonstrates:
- International cooperation in computational biology
- Multi-institutional validation of methodology
- Broad scientific impact across different research communities
- Standardization of computational approaches for NIPT data analysis

**Section sources**
- [README.md:17](file://README.md#L17)

## Methodology Acknowledgment Guidelines

### Proper Attribution Practices

When reporting results obtained through BaseVar2, researchers should acknowledge both the tool and the underlying methodology:

**Required acknowledgment elements:**
1. **Tool citation:** BaseVar2 software tool
2. **Methodology:** Ultra-low-depth WGS variant calling methodology
3. **Version information:** Specific BaseVar2 version used
4. **Computational approach:** Maximum likelihood and likelihood ratio models

### Research Paper Reporting Standards

**Recommended acknowledgment text:**
- "We acknowledge the use of BaseVar2 v2.0.0 (Huang et al., 2024) for variant calling from ultra-low-depth WGS data, employing maximum likelihood and likelihood ratio models as described in the original publication."

### Conference Presentation Acknowledgments

For conference presentations:
- Include software citation in slide footers
- Reference the specific methodology used
- Acknowledge computational resources if applicable
- Provide contact information for software inquiries

### Grant and Funding Acknowledgments

When acknowledging funding sources:
- Reference grant numbers supporting the research
- Acknowledge institutional support for computational infrastructure
- Credit funding agencies for supporting the development of BaseVar2

**Section sources**
- [README.md:9](file://README.md#L9)
- [src/basetype.cpp:170-171](file://src/basetype.cpp#L170-L171)

## Technical Implementation Notes

### Software Architecture Implications

The BaseVar2 implementation reflects the computational sophistication required for NIPT data analysis:

- **C++17 implementation:** Demonstrates high-performance computing requirements
- **Thread pool architecture:** Supports parallel processing for large-scale genomic analysis
- **Memory optimization:** Addresses computational constraints of ultra-low-depth data processing
- **Integration with HTSlib:** Leverages established bioinformatics infrastructure

### Algorithmic Foundations

The methodology employs sophisticated statistical approaches:
- **Maximum likelihood estimation:** For accurate variant detection
- **Likelihood ratio testing:** For statistical significance assessment
- **EM algorithm implementation:** For population allele frequency estimation
- **Quality score calculations:** For variant confidence assessment

**Section sources**
- [CMakeLists.txt:5](file://CMakeLists.txt#L5)
- [src/algorithm.h:174-177](file://src/algorithm.h#L174-L177)
- [src/variant_caller.h:33](file://src/variant_caller.h#L33)

## Conclusion

Proper academic attribution for BaseVar2 research requires recognition of both the software tool and the underlying methodology. The primary publication in Cell Genomics (2024) establishes the foundation for using NIPT sequencing data in human genetic investigation, while BaseVar2 provides the computational framework for implementing these methodologies.

Researchers utilizing BaseVar2 should follow the citation guidelines outlined above, ensuring proper acknowledgment of the collaborative research team, institutional affiliations, and the specific version of the software used. This approach maintains scientific integrity and supports continued development of computational tools for genomic research.

The combination of rigorous methodology, international collaboration, and open-source software development demonstrated by BaseVar2 represents a model for advancing computational genomics research while maintaining proper academic standards.