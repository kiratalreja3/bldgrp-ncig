#!/bin/bash

input_dir=/g/data/te53/sj2852/blgrp/analyses/flye_assemblies/NEW_FLYE/hapdup/output_files
ref_file=/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/NEW_COORDS/chm13/blgrp_coords_chm13_anchorflank.fa
outdir=/g/data/te53/sj2852/blgrp/analyses/minimap_flye_assemblies/NEW

for subdir in "$input_dir"/*/; do

    folder_name=$(basename "$subdir")

    phased_1="$subdir/hapdup_phased_1.fasta"
    phased_2="$subdir/hapdup_phased_2.fasta"

    if [[ -f "$phased_1" && -f "$phased_2" ]]; then

        qsub -v INPUT_FASTA="$phased_1",INPUT_REF="$ref_file",OUTPUT_DIR="$outdir",OUTPUT_FILE="${folder_name}_phased_1.paf" \
        /g/data/te53/sj2852/blgrp/scripts/run_and_submit_scripts/minimap_contigs_new.sh

        qsub -v INPUT_FASTA="$phased_2",INPUT_REF="$ref_file",OUTPUT_DIR="$outdir",OUTPUT_FILE="${folder_name}_phased_2.paf" \
        /g/data/te53/sj2852/blgrp/scripts/run_and_submit_scripts/minimap_contigs_new.sh

    else
        echo "Either hapdup_phased_1.fasta or hapdup_phased_2.fasta not found in $subdir"

    fi

done
