#!/bin/bash

# Print the first 4 lines (first read) of each fastq.gz file in the current directory

for file in *.fastq.gz; do
    if [ -f "$file" ]; then
        echo "========================================="
        echo "File: $file"
        echo "========================================="
        zcat "$file" | head -n 4
        echo ""
    fi
done

echo "Done!"

