#!/bin/bash
#PBS -l ncpus=8,mem=32GB,walltime=00:30:00,storage=gdata/te53+gdata/if89
#PBS -N bamextract
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/
#PBS -q normal
#PBS -P te53

#this script generates the files found in the assembly_stats folder which are used by the /g/data/te53/sj2852/blgrp/scripts/plotting/quality_read_depth_read_length_plotting.R script

module load samtools

bam_folder="/g/data/te53/sj2852/blgrp/data/KGPdata/onekgp100plusgrch38_bam"


output_file="/g/data/te53/sj2852/blgrp/metadata/assembly_stats/newer_stats/grch38_combined_assembly_stats_longform.tsv"

> $output_file

for bam_file in "$bam_folder"/*.bam; do
    [[ "$bam_file" == *.bam.bai ]] && continue

    sample_name=$(basename "$bam_file" .bam)

    samtools stats "$bam_file" | grep ^SN | cut -f 2- | awk -F ': ' -v sample="$sample_name" '{print sample "\t" $1 "\t" $2}' >> $output_file
done

#transpose the data in the output.tsv file into new columns

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

