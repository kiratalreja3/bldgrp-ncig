#!/bin/bash
#PBS -l ncpus=4,mem=16GB,walltime=04:00:00,storage=gdata/te53+gdata/if89
#PBS -N mappingq
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/
#PBS -q normalbw

#Author: Sarah Jackson 

# this code generates mapping quality statistics at a sample and region level

module load samtools

#For Sample Level

#BAM_DIR="/g/data/te53/sj2852/blgrp/analyses/reference_alignment_initial/NEW_alignments/NCIG_all/ncig_bam_and_bai" 

#OUTPUT_DIR="/g/data/te53/sj2852/blgrp/metadata/assembly_stats/NEW/mapping_quality/sample_level"

#for bam_file in "$BAM_DIR"/*.bam; do
    #[[ "$bam_file" == *.bam.bai ]] && continue

    #sample_name=$(basename "$bam_file" .bam)
    
    #file_name=$(basename "$bam_file" .bam)
    #output_file="$OUTPUT_DIR/${file_name}_mapping_quality.txt"

    #samtools view "$bam_file" | awk '{print $5}' | sort | uniq -c >"$output_file"

#done 

#For Region Level

#parent_directory="/g/data/te53/sj2852/blgrp/metadata/assembly_stats/NEW/by_blgrp_region/region_bam_files"

#output_base_folder="/g/data/te53/sj2852/blgrp/metadata/assembly_stats/NEW/mapping_quality/region_level"

#input_directories=$(find "$parent_directory" -mindepth 1 -maxdepth 1 -type d)

# Loop over each found subdirectory
#for bam_folder in $input_directories; do

    # Get the name of the current subdirectory
    #subdirectory=$(basename "$bam_folder")

    # Create an output folder for the current subdirectory
    #output_folder="${output_base_folder}/${subdirectory}"
    #mkdir -p "$output_folder"

    #Run on BAM files
   #for bam_file in "$bam_folder"/*.bam; do
        #[[ "$bam_file" == *.bam.bai ]] && continue

        #sample_name=$(basename "$bam_file" .bam)

        # Set the output file name based on the sample name
        #output_file="${output_folder}/${sample_name}_mapping_quality.txt"

        # Calculate mapping quality and save to the output file
        #samtools view "$bam_file" | awk '{print $5}' | sort | uniq -c > "$output_file"
    #done

#done

#Calculate Average Mapping Quality Per Region

reg_dir="/g/data/te53/sj2852/blgrp/metadata/assembly_stats/NEW/mapping_quality/region_level"

output_file="/g/data/te53/sj2852/blgrp/metadata/assembly_stats/NEW/mapping_quality/average_mapping_region.txt"

> "$output_file"

for dir in "$reg_dir"/*; do
    if [ -d "$dir" ]; then

     subdirectory=$(basename "$dir")

    total_quality=0
    total_count=0

     # Loop through each file in the subdirectory
    for file in "$dir"/*.txt; do
      # Read the mapping quality and count from each file
      while read count quality; do
        total_quality=$((total_quality + (count * quality)))
        total_count=$((total_count + count))
      done < "$file"
    done

    if [ $total_count -gt 0 ]; then
      average_quality=$(echo "$total_quality / $total_count" | bc -l)
    else
      average_quality=0
    fi

   # Append the result for this subdirectory to the output file
    echo "$subdirectory: $average_quality" >> "$output_file"
  fi
done