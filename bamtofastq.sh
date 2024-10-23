#!/bin/bash
#PBS -l ncpus=4,mem=16GB,walltime=04:00:00,storage=gdata/te53+gdata/if89
#PBS -N bamtofastq_chm13
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/
#PBS -q normalbw

input_folder="/g/data/te53/sj2852/blgrp/data/NEW_KGP_DATA/onekgpchm13_bam"
output_folder="/g/data/te53/sj2852/blgrp/data/NEW_KGP_DATA/kgp_data_all_fastq"

module load samtools 

# Iterate over each BAM file in the input folder
for input_bam in "${input_folder}"/*.bam; do
    # Check if the current item in the loop is a regular file
    if [ -f "$input_bam" ]; then
        # Define output file path for FASTQ file
        output_fastq="${output_folder}/$(basename "${input_bam}" .bam).fastq.gz"
        
        # Generate FASTQ file from BAM input and compress it
        echo "Converting BAM to FASTQ: ${input_bam} -> ${output_fastq}"
        samtools fastq "${input_bam}" | gzip > "${output_fastq}"
    fi
done