#!/bin/bash
#PBS -l ncpus=4,mem=16GB,walltime=08:00:00,storage=gdata/te53+gdata/if89
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/
#PBS -q normalbw
#PBS -N samtools_stats

#Author: Sarah Jackson & Kirat Alreja

#this code creates samtools stats and transposed files for each blood group region
#/g/data/te53/sj2852/blgrp/metadata/assembly_stats/NEW/by_blgrp_region/stats_files

# Load samtools module
module load samtools

parent_directory="/g/data/te53/sj2852/blgrp/metadata/assembly_stats/NEW/by_blgrp_region/region_bam_files"

output_folder="/g/data/te53/sj2852/blgrp/metadata/assembly_stats/NEW/by_blgrp_region/stats_files"

input_directories=$(find "$parent_directory" -mindepth 1 -maxdepth 1 -type d)

# Loop over each found subdirectory
for bam_folder in $input_directories; do

    # Get the name of the current subdirectory
    subdirectory=$(basename "$bam_folder")

    # Set the output file name based on the subdirectory name
    output_file="${output_folder}/${subdirectory}_stats.tsv"

    if [[ -f "$output_file" ]]; then
        echo "Output file for $subdirectory already exists. Skipping this directory."
        continue
    fi

    # Clear or create the output file
    > "$output_file"

    # Process each BAM file in the subdirectory
    for bam_file in "$bam_folder"/*.bam; do
        [[ "$bam_file" == *.bam.bai ]] && continue

        sample_name=$(basename "$bam_file" .bam)

        # Run samtools stats and append the results to the output file
        samtools stats "$bam_file" | grep ^SN | cut -f 2- | awk -F ': ' -v sample="$sample_name" '{print sample "\t" $1 "\t" $2}' >> "$output_file"
    done

done

# Loop to process each output file and generate a transposed version
for output_file in "$output_folder"/*_stats.tsv; do
    awk '
    BEGIN {
    FS = "\t"  # Set the field separator to tab
    OFS = "\t"  # Set the output field separator to tab
    header = "Sample" # Initialize the header variable with "Sample"
}

{
    # Store the sample ID (first column) in the sample variable
    sample = $1

    # Store the statistic name (second column) in the stat variable
    stat = $2

    # Store the statistic value (third column) in the value variable
    value = $3

    # Create an array to hold the data for each sample
    # The array key is the sample ID and the value is a string containing statistics
    data[sample] = data[sample] OFS value

    # Collect unique statistic names for the header
    if (!(stat in header_fields)) {
    header = header OFS stat
    header_fields[stat] = 1
    }
}

# End of each record, print the collected data for each sample
END {
    # Iterate over each sample ID in the data array
    for (sample in data) {
        # Print the sample ID followed by the collected statistics
        print sample, data[sample]
    }
}
    ' "$output_file" > "${output_file%.tsv}_transposed.tsv"
done