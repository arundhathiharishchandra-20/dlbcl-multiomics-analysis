# DLBCL Multiomics Analysis

## Overview
This repository contains scripts for a multi-omics analysis workflow of Diffuse Large B-Cell Lymphoma (DLBCL), integrating RNA sequencing and whole-genome bisulfite sequencing (WGBS) data via Multi-Omics Factor Analysis (MOFA2).

## Repository Structure
RNA_seq/    - RNA-seq preprocessing, WGCNA, pathway enrichment
WGBS/       - WGBS alignment, methylation quantification, gene mapping
multiomics/ - Correlation analysis + MOFA2 integration
reference/  - Genome preparation scripts

## Pipeline Order
1. RNA_seq/download.sh -> run_fastp_batch*.sh -> HISAT2 -> featureCounts
2. RNA_seq/wgcna_outputs.R -> WGCNA modules + hub genes
3. RNA_seq/pathway_analysis.R -> GO/KEGG enrichment
4. WGBS/download_wgbs.sh -> bismark_alignment.sh -> step3_deduplication.sh -> step4_methylation_extraction.sh
5. WGBS/methylkit_analysis.R -> tile-based methylation quantification
6. WGBS/create_gene_methylation_methylkit.py -> gene-level methylation matrix
7. multiomics/multiomics_updated.R -> correlation, GO enrichment, STRING network
8. multiomics/mofa_updated.R -> MOFA2 integration

## Software
R (>=4.3), Bioconductor (methylKit, MOFA2, edgeR, WGCNA, gprofiler2, STRINGdb, org.Hs.eg.db), Python 3, SRA Toolkit, Bismark, Bowtie2, Samtools, FastQC, fastp, HISAT2

## Status
RNA-seq, WGBS, and multi-omics integration (correlation + MOFA2) pipelines are complete. Part of an ongoing M.Sc. dissertation project.
