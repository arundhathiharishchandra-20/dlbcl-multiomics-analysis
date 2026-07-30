#!/bin/bash
BASE_DIR="/mnt/d/lymphoma_project"   # <-- change this to your own path
set -euo pipefail

OUTDIR="$BASE_DIR/WGBS/merged_raw"
TRIMDIR="$BASE_DIR/WGBS/merged_trimmed"
ALIGNDIR="$BASE_DIR/WGBS/merged_aligned"
METHYLDIR="$BASE_DIR/WGBS/merged_methylation"
LOGDIR="$BASE_DIR/WGBS/logs"
TMPDIR="$BASE_DIR/WGBS/tmp"
GENOME="$BASE_DIR/reference/genome_fasta"
BISMARK="$BASE_DIR/Bismark/bismark"
METHEXT="$BASE_DIR/Bismark/bismark_methylation_extractor"
THREADS=2
SAMPLE_MAP="$BASE_DIR/WGBS/gsm_51_paired_samples.txt"

mkdir -p "$OUTDIR" "$TRIMDIR" "$ALIGNDIR" "$METHYLDIR" "$LOGDIR" "$TMPDIR"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOGDIR/pipeline.log"; }

# Check if sample is already successfully processed
is_successful() {
local sample="$1"
local base="${METHYLDIR}/${sample}_merged.deduplicated"
[[ -s "${base}.bismark.cov.gz" && -s "${base}.CpG_report.txt.gz" ]] || return 1
local cpg
cpg=$(zcat "${base}.bismark.cov.gz" 2>/dev/null | wc -l || echo 0)
(( cpg > 500000 ))
}

while IFS=$'\t' read -r GSM FIRST_RUN EXTRA_RUNS || [[ -n "$GSM" ]]; do
[[ -z "$GSM" || "$GSM" =~ ^# ]] && continue

if is_successful "$GSM"; then
log "Skipping $GSM (already complete)"
continue
fi

log "=== Starting sample: $GSM ==="
SAMPLE_DIR="$OUTDIR/$GSM"
mkdir -p "$SAMPLE_DIR"

if [[ -n "$EXTRA_RUNS" ]]; then
ALL_SRRS="${FIRST_RUN},${EXTRA_RUNS}"
else
ALL_SRRS="${FIRST_RUN}"
fi

IFS=',' read -ra SRRS <<< "$ALL_SRRS"

# ── Step 1: Download, trim, align each SRR separately ─────────────
BAM_LIST=()

for SRR in "${SRRS[@]}"; do
SRR="${SRR// /}"
[[ -z "$SRR" ]] && continue

SRR_BAM="$ALIGNDIR/${GSM}_${SRR}_bismark_bt2.bam"

# Skip if this SRR BAM already done
if [[ -s "$SRR_BAM" ]]; then
log "$SRR BAM already exists - skipping"
BAM_LIST+=("$SRR_BAM")
continue
fi

# Download if fastq missing
if [[ ! -s "$SAMPLE_DIR/${SRR}.fastq" ]]; then
rm -f "$SAMPLE_DIR/${SRR}/${SRR}.sra.lock" "$SAMPLE_DIR/${SRR}/${SRR}.sra.tmp" "$SAMPLE_DIR/${SRR}/${SRR}.sra.prf"
log "Downloading $SRR ..."
prefetch "$SRR" --output-directory "$SAMPLE_DIR" 2>&1 | tee -a "$LOGDIR/${GSM}.log" || { log "WARNING: prefetch failed for $SRR - skipping"; continue; }
SRA_FILE="$SAMPLE_DIR/${SRR}/${SRR}.sra"
if [[ -s "$SRA_FILE" ]]; then
log "Converting $SRR to fastq..."
fasterq-dump "$SRA_FILE" --outdir "$SAMPLE_DIR" --temp "$TMPDIR" --split-files --threads 4 2>&1 | tee -a "$LOGDIR/${GSM}.log" || { log "WARNING: fasterq-dump failed for $SRR - skipping"; continue; }
if [[ -s "$SAMPLE_DIR/${SRR}_1.fastq" && ! -s "$SAMPLE_DIR/${SRR}.fastq" ]]; then
mv "$SAMPLE_DIR/${SRR}_1.fastq" "$SAMPLE_DIR/${SRR}.fastq"
rm -f "$SAMPLE_DIR/${SRR}_2.fastq"
fi
rm -rf "$SAMPLE_DIR/${SRR}"
else
log "WARNING: SRA file not found for $SRR - skipping"
continue
fi
fi

[[ -s "$SAMPLE_DIR/${SRR}.fastq" ]] || { log "WARNING: No fastq for $SRR - skipping"; continue; }

# Trim
SRR_TRIMMED="$TRIMDIR/${GSM}_${SRR}_trimmed.fastq"
if [[ ! -s "$SRR_TRIMMED" ]]; then
log "Trimming $SRR ..."
fastp -i "$SAMPLE_DIR/${SRR}.fastq" -o "$SRR_TRIMMED" --thread "$THREADS" --qualified_quality_phred 20 --length_required 36 2>&1 | tee -a "$LOGDIR/${GSM}.log"
fi
[[ -s "$SRR_TRIMMED" ]] || { log "WARNING: fastp failed for $SRR - skipping"; continue; }

# Align
log "Aligning $SRR ..."
"$BISMARK" --genome "$GENOME" --output_dir "$ALIGNDIR" -p 2 "$SRR_TRIMMED" 2>&1 | tee -a "$LOGDIR/${GSM}.log"

# Bismark names output as: inputfilename_bismark_bt2.bam
RAW_BAM="$ALIGNDIR/${GSM}_${SRR}_trimmed_bismark_bt2.bam"
if [[ -s "$RAW_BAM" ]]; then
mv "$RAW_BAM" "$SRR_BAM"
BAM_LIST+=("$SRR_BAM")
log "$SRR aligned successfully"
else
log "WARNING: Bismark BAM missing for $SRR - skipping"
fi

# Cleanup trimmed fastq and raw fastq immediately
rm -f "$SRR_TRIMMED"
rm -f "$SAMPLE_DIR/${SRR}.fastq"
rm -f "$TMPDIR"/*

done

# Check we have BAMs to merge
if (( ${#BAM_LIST[@]} == 0 )); then
log "WARNING: No BAM files for $GSM - skipping sample"
continue
fi
log "Found ${#BAM_LIST[@]} BAM files for $GSM - merging..."

# ── Step 2: Merge all SRR BAMs ────────────────────────────────────
MERGED_BAM="$ALIGNDIR/${GSM}_merged.bam"
if (( ${#BAM_LIST[@]} == 1 )); then
cp "${BAM_LIST[0]}" "$MERGED_BAM"
else
samtools merge -f "$MERGED_BAM" "${BAM_LIST[@]}" 2>&1 | tee -a "$LOGDIR/${GSM}.log"
fi
[[ -s "$MERGED_BAM" ]] || { log "WARNING: BAM merge failed for $GSM - skipping"; continue; }

# Delete individual SRR BAMs
for bam in "${BAM_LIST[@]}"; do
rm -f "$bam"
done
log "BAMs merged for $GSM"

# ── Step 3: Deduplicate with samtools markdup ──────────────────────
log "Deduplicating $GSM ..."
samtools sort -n "$MERGED_BAM" -o "${MERGED_BAM%.bam}_nsorted.bam" 2>&1 | tee -a "$LOGDIR/${GSM}.log"
samtools fixmate -m "${MERGED_BAM%.bam}_nsorted.bam" "${MERGED_BAM%.bam}_fixmate.bam" 2>&1 | tee -a "$LOGDIR/${GSM}.log"
samtools sort "${MERGED_BAM%.bam}_fixmate.bam" -o "${MERGED_BAM%.bam}_sorted.bam" 2>&1 | tee -a "$LOGDIR/${GSM}.log"
samtools markdup -r "${MERGED_BAM%.bam}_sorted.bam" "${MERGED_BAM%.bam}.deduplicated.bam" 2>&1 | tee -a "$LOGDIR/${GSM}.log"
rm -f "${MERGED_BAM%.bam}_nsorted.bam" "${MERGED_BAM%.bam}_fixmate.bam" "${MERGED_BAM%.bam}_sorted.bam"
DEDUP_BAM="${MERGED_BAM%.bam}.deduplicated.bam"
[[ -s "$DEDUP_BAM" ]] || { log "WARNING: Deduplication failed for $GSM - skipping"; continue; }
rm -f "$MERGED_BAM"
log "Deduplication done for $GSM"

# ── Step 4: Methylation extraction ────────────────────────────────
log "Extracting methylation for $GSM ..."
"$METHEXT" --single-end --comprehensive --bedGraph --counts --gzip --cytosine_report --genome_folder "$GENOME" --parallel 2 --output "$METHYLDIR" "$DEDUP_BAM" 2>&1 | tee -a "$LOGDIR/${GSM}.log"

# Validate output
# Validate output
BASE="${METHYLDIR}/${GSM}_merged.deduplicated"
[[ -s "${BASE}.bismark.cov.gz" ]] || { log "WARNING: Coverage file missing for $GSM - skipping"; continue; }
[[ -s "${BASE}.CpG_report.txt.gz" ]] || { log "WARNING: CpG report missing for $GSM - skipping"; continue; }

# Final cleanup
rm -f "$DEDUP_BAM"
rm -f "$TMPDIR"/*
log "$GSM finished successfully!"

done < "$SAMPLE_MAP"

log "=== Pipeline completed successfully! ==="
