#!/bin/bash
BASE_DIR="/mnt/d/lymphoma_project"   # <-- change this to your own path

FASTQ_DIR="$BASE_DIR/raw_fastq"
OUT_DIR="/mnt/c/Users/Admin/Downloads/lymphoma_project/fastp_output"
LIST="$BASE_DIR/downloaded_list.txt"
THREADS=4
COUNT=81

sed -n '81,95p' "$LIST" | while read SAMPLE; do
    NUM=$(printf "%03d" $COUNT)
    echo ">>> Processing [$NUM]: $SAMPLE"

    fastp \
        -i "$FASTQ_DIR/${SAMPLE}_1.fastq" \
        -I "$FASTQ_DIR/${SAMPLE}_2.fastq" \
        -o "$OUT_DIR/clean_reads/TRIM${NUM}_1.fastq.gz" \
        -O "$OUT_DIR/clean_reads/TRIM${NUM}_2.fastq.gz" \
        -h "$OUT_DIR/reports/${NUM}.html" \
        -j "$OUT_DIR/reports/${NUM}.json" \
        --thread $THREADS \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 36 \
        --correction \
        2>> "$OUT_DIR/fastp_run.log"

    echo "<<< Done: $SAMPLE → TRIM${NUM}"
    COUNT=$((COUNT + 1))
done

echo "BATCH 3 DONE! ALL 95 SAMPLES COMPLETE!"
