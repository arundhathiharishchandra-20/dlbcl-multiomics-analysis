#!/bin/bash
BASE_DIR="/mnt/d/lymphoma_project"   # <-- change this to your own path
METHYLDIR="$BASE_DIR/WGBS/merged_methylation"
SAMPLE_MAP="$BASE_DIR/WGBS/gsm_srr_mapping_complete.txt"
OUTFILE="$BASE_DIR/WGBS/good_samples.txt"

> "$OUTFILE"

while IFS=$'\t' read -r GSM FIRST_RUN EXTRA_RUNS || [[ -n "$GSM" ]]; do
    [[ -z "$GSM" || "$GSM" =~ ^# ]] && continue
    base="${METHYLDIR}/${GSM}_trimmed_bismark_bt2.deduplicated"
    [[ ! -s "${base}.cov.gz" || ! -s "${base}.CpG_report.txt.gz" ]] && continue
    cpg=$(zcat "${base}.cov.gz" 2>/dev/null | wc -l || echo 0)
    if (( cpg > 500000 )); then
        echo "$GSM" >> "$OUTFILE"
    fi
done < "$SAMPLE_MAP"

echo "Good samples saved to $OUTFILE"
cat "$OUTFILE"
