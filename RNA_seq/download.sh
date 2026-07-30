#!/bin/bash
BASE_DIR="/mnt/d/lymphoma_project"   # <-- change this to your own path
cd $BASE_DIR/raw_fastq

while read srr; do
    if [ ! -f "${srr}_1.fastq" ]; then
        echo "Downloading $srr..."
        prefetch $srr -O ./
        fasterq-dump $srr --split-files -p \
            -O $BASE_DIR/raw_fastq/
        rm -rf $BASE_DIR/raw_fastq/${srr}
        echo "Done: $srr"
    else
        echo "Skipping $srr - already downloaded!"
    fi
done < "$(dirname "$0")/accessions.txt"
