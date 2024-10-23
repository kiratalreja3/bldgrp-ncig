#!/bin/bash
#PBS -l ncpus=4,mem=16GB,walltime=04:00:00,storage=gdata/te53+gdata/if89
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/
#PBS -q normalbw
#PBS -N samtools_stats

#Author: Sarah Jackson & Kirat Alreja

#this code creates samtools stats and transposed files per sample

# Load samtools module
module load samtools

bam_folder="/g/data/te53/sj2852/blgrp/analyses/reference_alignment_initial/NEW_alignments/NCIG_all/ncig_bam_and_bai"

# Set the output folder
output_folder="/g/data/te53/sj2852/blgrp/metadata/assembly_stats/NEW"

# Make sure the output folder exists
mkdir -p "$output_folder"

# Set the output file name
output_file="${output_folder}/NCIG_combined_assembly_stats_longform.tsv"

# Clear or create the output file
> "$output_file"

# Process each BAM file in the specified directory
for bam_file in "$bam_folder"/*.bam; do
    # Skip BAM index files
    [[ "$bam_file" == *.bam.bai ]] && continue

    sample_name=$(basename "$bam_file" .bam)

    # Run samtools stats and append the results to the output file
    samtools stats "$bam_file" | grep ^SN | cut -f 2- | awk -F ': ' -v sample="$sample_name" '{print sample "\t" $1 "\t" $2}' >> "$output_file"
done

# Process the output file to generate a transposed version
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
    # Print the header
    print header
    # Iterate over each sample ID in the data array
    for (sample in data) {
        # Print the sample ID followed by the collected statistics
        print sample, data[sample]
    }
}
' "$output_file" > "${output_file%.tsv}_transposed.tsv"