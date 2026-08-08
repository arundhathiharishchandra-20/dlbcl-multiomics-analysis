"""
Create Gene_Methylation.csv from methylKit CpG matrix
Input:  CpG_matrix_for_gene_meth.csv (from methylkit_analysis.R)
Output: Gene_Methylation.csv (for multiomics.R and mofa.R)
"""

import gzip
import pandas as pd
import numpy as np
from collections import defaultdict

METHYLKIT_OUT = "/mnt/d/lymphoma_project/methylkit_output"
REFGENE       = "/mnt/d/lymphoma_project/WGBS/refGene.txt.gz"
OUTFILE       = "/mnt/d/lymphoma_project/RNA_seq/counts/Gene_Methylation.csv"
CpG_MATRIX    = f"{METHYLKIT_OUT}/CpG_matrix_for_gene_meth.csv"

print("Step 1: Loading CpG methylation matrix...")
cpg_df = pd.read_csv(CpG_MATRIX, index_col=0)
print(f"CpG matrix: {cpg_df.shape[0]} CpGs x {cpg_df.shape[1]} samples")

# Parse chr:pos from index
print("Step 2: Parsing CpG coordinates...")
chroms = []
positions = []
for idx in cpg_df.index:
    parts = idx.split(":")
    chroms.append(parts[0])
    positions.append(int(parts[1]))

cpg_df['chrom'] = chroms
cpg_df['pos'] = positions

print("Step 3: Loading gene promoter coordinates...")
promoters = defaultdict(list)
with gzip.open(REFGENE, 'rt') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 13:
            continue
        chrom     = parts[2]
        strand    = parts[3]
        tx_start  = int(parts[4])
        tx_end    = int(parts[5])
        gene_name = parts[12]

        if strand == '+':
            prom_start = max(0, tx_start - 2000)
            prom_end   = tx_start + 500
        else:
            prom_start = max(0, tx_end - 500)
            prom_end   = tx_end + 2000

        promoters[chrom].append((prom_start, prom_end, gene_name))

print(f"Loaded promoters for {len(promoters)} chromosomes")

print("Step 4: Overlapping CpGs with promoters...")
sample_cols = [c for c in cpg_df.columns if c not in ['chrom', 'pos']]
gene_meth   = defaultdict(lambda: {s: [] for s in sample_cols})
matched     = 0

for idx, row in cpg_df.iterrows():
    chrom = row['chrom']
    pos   = row['pos']

    # Add chr prefix if needed
    chrom_key = "chr" + chrom if not chrom.startswith("chr") else chrom

    if chrom_key not in promoters:
        continue

    for (prom_start, prom_end, gene_name) in promoters[chrom_key]:
        if prom_start <= pos <= prom_end:
            for s in sample_cols:
                val = row[s]
                if not pd.isna(val):
                    gene_meth[gene_name][s].append(val)
            matched += 1
            break

print(f"CpGs matched to promoters: {matched}")
print(f"Genes with methylation data: {len(gene_meth)}")

print("Step 5: Averaging methylation per gene...")
result = {}
for gene, sample_vals in gene_meth.items():
    result[gene] = {}
    for s in sample_cols:
        vals = sample_vals[s]
        result[gene][s] = np.mean(vals) if vals else np.nan

df = pd.DataFrame(result).T
df.index.name = 'Gene'

# Remove genes with too many NAs
na_frac = df.isna().mean(axis=1)
df = df[na_frac < 0.5]
print(f"Genes after NA filter (>50% NA removed): {df.shape[0]}")

print("Step 6: Saving Gene_Methylation.csv...")
df.to_csv(OUTFILE)
print(f"Saved to {OUTFILE}")
print(f"Final matrix: {df.shape[0]} genes x {df.shape[1]} samples")
print("Done! Now run multiomics_updated.R and mofa_updated.R")
