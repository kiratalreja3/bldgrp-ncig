#!/bin/bash
#PBS -l ncpus=8,mem=32GB,walltime=10:00:00,storage=gdata/te53+gdata/if89
#PBS -N minimap_contigs
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/
#PBS -q normal
#PBS -P te53

module load minimap2/2.28

INPUT_FASTA=$INPUT_FASTA
INPUT_REF=$INPUT_REF
OUTPUT_DIR=$OUTPUT_DIR
OUTPUT_FILE=$OUTPUT_FILE

minimap2 --cs --secondary=no "$INPUT_REF" "$INPUT_FASTA" > "$OUTPUT_DIR/$OUTPUT_FILE"

#Other way round
#minimap2 --cs --secondary=no "$INPUT_FASTA" "$INPUT_REF" > "$OUTPUT_DIR/$OUTPUT_FILE"

echo "Mapping complete. Output saved to $OUTPUT_FILE"

#grep 'tp:A:P' "$OUTPUT_DIR/$OUTPUT_FILE" | cut -f1,2,6,7,11,12 > "$OUTPUT_DIR/${OUTPUT_FILE%.paf}_filtered.paf"

#rm "$OUTPUT_DIR/$OUTPUT_FILE"
