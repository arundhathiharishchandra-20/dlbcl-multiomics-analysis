# DLBCL Multiomics Analysis

## Overview

This repository contains scripts for a multiomics analysis workflow of Diffuse Large B-Cell Lymphoma (DLBCL), integrating RNA sequencing and whole-genome bisulfite sequencing (WGBS) data.

## Repository Structure

```
RNA_seq/
    RNA-seq preprocessing and downstream analysis

WGBS/
    WGBS preprocessing, alignment, deduplication and methylation analysis

reference/
    Genome preparation scripts

mofa_updated.R
    Multi-Omics Factor Analysis (MOFA) workflow

multiomics_updated.R
    Integration of RNA-seq and DNA methylation analyses
```

## Software

- R (>= 4.3)
- Bismark
- Bowtie2
- Samtools
- FastQC
- fastp

## Status

This is an active research project. Some scripts (especially MOFA integration) are still under development and may be updated as the analysis progresses.# dlbcl-multiomics-analysis
RNA-seq and WGBS pipeline for multi-omics epigenetic analysis of DLBCL lymphoma
