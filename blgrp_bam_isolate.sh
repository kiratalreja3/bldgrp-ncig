#!/bin/bash
#PBS -l ncpus=4,mem=16GB,walltime=04:00:00,storage=gdata/te53+gdata/if89
#PBS -N bamextract_ncig
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/
#PBS -q normalbw

#Author: Sarah Jackson & Kirat Alreja 

#bam for blood groups
#this script creates the output files of bam reads per sample per blgrp region
#(/g/data/te53/sj2852/blgrp/metadata/assembly_stats/by_blgrp_region/region_bam_files)

module load samtools


#step 1 - extract bam for each blgrp region for each sample
BED_FILE="/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/NEW_COORDS/chm13/blgrp_coords_chm13_anchorflank_named.bed"
BAM_FILES_PATH="/g/data/te53/sj2852/blgrp/analyses/reference_alignment_initial/NEW_alignments/NCIG_all/ncig_bam_and_bai"
OUTPUT_BASE_DIR="/g/data/te53/sj2852/blgrp/metadata/assembly_stats/NEW/by_blgrp_region/region_bam_files"

while IFS=$'\t' read -r chr start end gene; do
    # Create the region and region_name variables
    end_adjusted=$((end + 1))
    region="${chr}:${start}-${end_adjusted}"
    region_name="${gene}"

    #check region
    echo "Processing region: $region" 

#Create the directory for the current region
    region_dir="${OUTPUT_BASE_DIR}/${region_name}"
    #mkdir -p "$region_dir"

 if [[ ! -d "$region_dir" ]]; then
        echo "Error: Output directory $region_dir does not exist. Skipping region $region_name."
        continue
    fi

for bam_file in "${BAM_FILES_PATH}"/*.bam; do

    #does BAM exist?
    if [[ ! -f "$bam_file" ]]; then
            echo "Warning: BAM file '$bam_file' does not exist."
            continue
        fi

    # Extract the sample name (assuming the sample name is the file name without the .bam extension)
    sample_name=$(basename "$bam_file" .bam)

    # Define the output BAM file path
    output_bam="${region_dir}/${sample_name}.bam"

      if [[ -f "$output_bam" ]]; then
            echo "Skipping existing file: $output_bam"
            continue
        fi

    if ! samtools view -b "$bam_file" "$region" > "$output_bam"; then
            echo "Error: Failed to process BAM file '$bam_file' for region '$region'."
    else
            echo "Successfully created: $output_bam"
    fi

   done
done < "$BED_FILE"


