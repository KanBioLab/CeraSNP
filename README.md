# CeraSNP: Two-Tier Core SNP Discovery Pipeline

[![Lab: KanBioLab](https://img.shields.io/badge/Laboratory-KanBioLab-084594.svg)](https://github.com/KanBioLab)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![R: >= 4.0](https://img.shields.io/badge/R-%3E%3D4.0-green.svg)](https://www.r-project.org/)

**Developed by Kan Laboratory (KanBioLab)**  
*Plant Genetics, Organelle Genomics & Computational Biology*

---

> **Associated Publication:**  
> *Comparative Plastomics Reveals Maternal Lineage Structure and SNP-Based Fingerprints of Ornamental Cherry Cultivars*

## Overview
**CeraSNP** is an end-to-end, two-tier cascaded adaptive feature selection pipeline designed for deterministic cultivar identification and germplasm fingerprinting:
- **Tier 1 (Plastome Backbone):** Forward adaptive greedy selection driven by marginal haplotype gain saturation.
- **Tier 2 (Nuclear 45S rDNA):** Residual collision minimization conditioned on maternal plastid backgrounds.
- **Assay QC & Validation:** In-silico KASP primer suitability assays coupled with strict zero-leakage Leave-One-Out Cross-Validation (LOOCV).

## Repository Structure
- `01_config_and_utils.R`: Global environment configuration, parameters, and core utilities.
- `02_greedy_engine.R`: Two-tier adaptive optimization engines and in-silico KASP QC functions.
- `03_discovery_execution.R`: Cohort discovery execution, marker selection, and comprehensive database construction.
- `04_loocv_validation.R`: Strict zero-leakage cross-validation and marker set retention/stability metrics.
- `05_visualization.R`: High-resolution vector visualizations (Trajectories, Hamming distance heatmap) and genomic polymorphism summaries.
- `run_pipeline.R`: One-click master controller.

## Dependencies
CeraSNP requires **R (>= 4.0)** and the following packages:
```R
# CRAN Packages
install.packages("pheatmap")

# Bioconductor Packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")
