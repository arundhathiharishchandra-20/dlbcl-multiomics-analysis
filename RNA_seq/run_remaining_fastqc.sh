#!/bin/bash
BASE_DIR="/mnt/d/lymphoma_project"   # <-- change this to your own path

INPUT_DIR="$BASE_DIR/raw_fastq"
OUTPUT_DIR="/mnt/c/Users/Admin/Downloads/lymphoma_project/fastqc_before"

mkdir -p "$OUTPUT_DIR"

for file in "$INPUT_DIR"/*.fastq
do
    base=$(basename "$file" .fastq)

    if [ -f "$OUTPUT_DIR/${base}_fastqc.html" ]; then
        echo "Skipping $base (already completed)"
    else
        echo "Running FastQC on $base"
        fastqc "$file" -o "$OUTPUT_DIR" -t 4
    fi
done

echo "All remaining FastQC analyses completed!"
