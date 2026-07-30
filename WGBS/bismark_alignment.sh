#!/bin/bash
BASE_DIR="/mnt/d/lymphoma_project"   # <-- change this to your own path

TRIM_DIR="$BASE_DIR/WGBS/fastp_output/clean_reads"
GENOME="$BASE_DIR/reference/genome_fasta"
BISMARK="$BASE_DIR/Bismark/bismark"
OUT_DIR="$BASE_DIR/WGBS/bismark_output/aligned"

mkdir -p $OUT_DIR

for SAMPLE in $TRIM_DIR/*.fastq.gz; do
    NAME=$(basename "$SAMPLE" ".fastq.gz")
    
    # Skip if already aligned
    if [ -f "$OUT_DIR/${NAME}_bismark_bt2.bam" ]; then
        echo "Already aligned, skipping: $NAME"
        continue
    fi

    echo ">>> Aligning: $NAME"

    $BISMARK \
    --genome $GENOME \
    --parallel 1 \
    --bowtie2 \
    -p 4 \
    -o $OUT_DIR \
    $SAMPLE

    echo "<<< Done: $NAME"
done

echo "ALL SAMPLES ALIGNED!"
