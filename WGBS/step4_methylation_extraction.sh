#!/bin/bash
BASE_DIR="/mnt/d/lymphoma_project"   # <-- change this to your own path

DEDUP_DIR="$BASE_DIR/WGBS/bismark_output/deduplicated"
METHYL_DIR="$BASE_DIR/WGBS/bismark_output/methylation"
GENOME="$BASE_DIR/reference/genome_fasta"
BISMARK_DIR="$BASE_DIR/Bismark"

mkdir -p $METHYL_DIR

echo "=============================="
echo "STEP 4: Methylation Extraction"
echo "=============================="

for BAM in $DEDUP_DIR/*.deduplicated.bam; do
    NAME=$(basename "$BAM" ".deduplicated.bam")

    if [ -f "$METHYL_DIR/${NAME}.bismark.cov.gz" ]; then
        echo "Already extracted, skipping: $NAME"
        continue
    fi

    echo ">>> Extracting methylation: $NAME"

    $BISMARK_DIR/bismark_methylation_extractor \
        --single-end \
        --comprehensive \
        --CX_context \
        --parallel 4 \
        --genome_folder $GENOME \
        --output $METHYL_DIR \
        --bedGraph \
        --counts \
        $BAM

    echo "<<< Done: $NAME"
done

echo "METHYLATION EXTRACTION COMPLETE!"
echo "CpG files saved in: $METHYL_DIR"
