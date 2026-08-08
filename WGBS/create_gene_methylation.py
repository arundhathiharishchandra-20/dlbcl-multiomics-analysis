import gzip
import os
import pandas as pd
import numpy as np
from collections import defaultdict

METHYLDIR = "/mnt/d/lymphoma_project/WGBS/merged_methylation"
REFGENE = "/mnt/d/lymphoma_project/WGBS/refGene.txt.gz"
OUTFILE = "/mnt/d/lymphoma_project/RNA_seq/counts/Gene_Methylation.csv"

print("Step 1: Loading gene promoter coordinates...")
promoters = {}
with gzip.open(REFGENE, 'rt') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 13:
            continue
        chrom = parts[2]
        strand = parts[3]
        tx_start = int(parts[4])
        tx_end = int(parts[5])
        gene_name = parts[12]
        if strand == '+':
            prom_start = max(0, tx_start - 2000)
            prom_end = tx_start + 500
        else:
            prom_start = max(0, tx_end - 500)
            prom_end = tx_end + 2000
        if chrom not in promoters:
            promoters[chrom] = []
        promoters[chrom].append((prom_start, prom_end, gene_name))

print(f"Loaded promoters for {len(promoters)} chromosomes")

print("Step 2: Finding coverage files...")
cov_files = sorted([f for f in os.listdir(METHYLDIR) 
                    if f.endswith('.bismark.cov.gz')])
sample_ids = [f.replace('_merged.deduplicated.bismark.cov.gz', '') 
              for f in cov_files]
print(f"Found {len(cov_files)} samples: {sample_ids}")

print("Step 3: Computing gene methylation per sample...")
gene_meth = defaultdict(lambda: defaultdict(list))

for i, (cov_file, sample_id) in enumerate(zip(cov_files, sample_ids)):
    print(f"  Processing {sample_id} ({i+1}/{len(cov_files)})...")
    filepath = os.path.join(METHYLDIR, cov_file)
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            chrom = "chr" + parts[0] if not parts[0].startswith("chr") else parts[0]
            pos = int(parts[1])
            try:
                methylated = int(parts[4])
                unmethylated = int(parts[5])
                total = methylated + unmethylated
                if total < 5:
                    continue
                pct_meth = methylated / total * 100
            except:
                continue
            if chrom not in promoters:
                continue
            for (prom_start, prom_end, gene_name) in promoters[chrom]:
                if prom_start <= pos <= prom_end:
                    gene_meth[gene_name][sample_id].append(pct_meth)

print(f"Step 4: Averaging methylation per gene...")
genes = sorted(gene_meth.keys())
result = {}
for gene in genes:
    result[gene] = {}
    for sample_id in sample_ids:
        values = gene_meth[gene].get(sample_id, [])
        result[gene][sample_id] = np.mean(values) if values else np.nan

df = pd.DataFrame(result).T
df.index.name = 'Gene'

print(f"Gene methylation matrix: {df.shape[0]} genes x {df.shape[1]} samples")

print("Step 5: Saving...")
df.to_csv(OUTFILE)
print(f"Saved to {OUTFILE}")
print("Done!")
