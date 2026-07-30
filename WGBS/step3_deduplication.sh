#!/bin/bash
BASE_DIR="/mnt/d/lymphoma_project"   # <-- change this to your own path

ALIGNED_DIR="$BASE_DIR/WGBS/bismark_output/aligned"
DEDUP_DIR="$BASE_DIR/WGBS/bismark_output/deduplicated"
BISMARK_DIR="$BASE_DIR/Bismark"

mkdir -p $DEDUP_DIR

echo "=============================="
echo "STEP 3: Deduplication"
echo "=============================="

for BAM in $ALIGNED_DIR/*_bismark_bt2.bam; do
    NAME=$(basename "$BAM" "_bismark_bt2.bam")

    if [ -f "$DEDUP_DIR/${NAME}_bismark_bt2.deduplicated.bam" ]; then
        echo "Already deduplicated, skipping: $NAME"
        continue
    fi

    echo ">>> Deduplicating: $NAME"

    $BISMARK_DIR/deduplicate_bismark \
        --single \
        --bam $BAM \
        --output_dir $DEDUP_DIR

    echo "<<< Done: $NAME"
done

echo "DEDUPLICATION COMPLETE!"
