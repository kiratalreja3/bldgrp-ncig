#!/bin/bash
#PBS -l ncpus=8,mem=32GB,walltime=00:30:00,storage=gdata/te53+gdata/if89
#PBS -N unzip_json
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/
#PBS -q normal
#PBS -P te53

# Directory containing the .json.gz files
input_dir="/g/data/te53/sj2852/blgrp/analyses/locityper_kgp/NEW/json_KGP"

# Loop over each .json.gz file in the directory
for gz_file in "$input_dir"/*.json.gz; do
    # Unzip the file
    gunzip "$gz_file"
done