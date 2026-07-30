#!/bin/bash
BASE_DIR="/mnt/d/lymphoma_project"   # <-- change this to your own path
METHYLDIR="$BASE_DIR/WGBS/merged_methylation"
SAMPLE_MAP="$BASE_DIR/WGBS/gsm_srr_mapping_complete.txt"

echo "=== FAILED SAMPLES ===" 
while IFS=$'\t' read -r GSM FIRST_RUN EXTRA_RUNS || [[ -n "$GSM" ]]; do
    [[ -z "$GSM" || "$GSM" =~ ^# ]] && continue
    base="${METHYLDIR}/${GSM}_trimmed_bismark_bt2.deduplicated"
    if [[ ! -s "${base}.cov.gz" || ! -s "${base}.CpG_report.txt.gz" ]]; then
        echo "MISSING FILES: $GSM"
        continue
    fi
    cpg=$(zcat "${base}.cov.gz" 2>/dev/null | wc -l || echo 0)
    if (( cpg <= 500000 )); then
        echo "LOW CpG ($cpg sites): $GSM"
    fi
done < "$SAMPLE_MAP"

echo ""
echo "=== PASSED SAMPLES ==="
while IFS=$'\t' read -r GSM FIRST_RUN EXTRA_RUNS || [[ -n "$GSM" ]]; do
    [[ -z "$GSM" || "$GSM" =~ ^# ]] && continue
    base="${METHYLDIR}/${GSM}_trimmed_bismark_bt2.deduplicated"
    [[ ! -s "${base}.cov.gz" ]] && continue
    cpg=$(zcat "${base}.cov.gz" 2>/dev/null | wc -l || echo 0)
    if (( cpg > 500000 )); then
        echo "OK ($cpg sites): $GSM"
    fi
done < "$SAMPLE_MAP"
