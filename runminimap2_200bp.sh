#!/bin/bash
#PBS -l ncpus=8,mem=32GB,walltime=03:00:00,storage=gdata/te53+gdata/if89
#PBS -N minimap_200bp
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/
#PBS -q normal
#PBS -P te53


set -e
set -o pipefail
set -u


module load minimap2/2.28
module load samtools

root_dir="/g/data/te53/sj2852/blgrp/analyses/flye_assemblies/hapdup"
output_dir="/g/data/te53/sj2852/blgrp/analyses/200bp_investigation/200bp_alignment_assemblies"
filtered_output_dir="/g/data/te53/sj2852/blgrp/analyses/200bp_investigation/200bp_alignment_filtered" 


for sample_dir in "$root_dir"/*/; do
  sample_id=$(basename "$sample_dir")

file1="${sample_dir}hapdup_phased_1.fasta"
file2="${sample_dir}hapdup_phased_2.fasta"

output1="${output_dir}/${sample_id}_hapdup_phased_1.minimap.paf"
output2="${output_dir}/${sample_id}_hapdup_phased_2.minimap.paf"


if [[ ! -f "$output1" ]]; then
    echo "Running minimap2 for $file1..."
    minimap2 --cs -t ${PBS_NCPUS} "$file1" "/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/chm13/200bp_anchors/bloodcoords_200bpflanks_chm13.fa" > "$output1"
  else
    echo "Minimap2 output for $file1 already exists. Skipping."
  fi

  if [[ ! -f "$output2" ]]; then
    echo "Running minimap2 for $file2..."
    minimap2 --cs -t ${PBS_NCPUS} "$file2" "/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/chm13/200bp_anchors/bloodcoords_200bpflanks_chm13.fa" > "$output2"
  else
    echo "Minimap2 output for $file2 already exists. Skipping."
  fi

echo "Minimap2 processing complete. Outputs saved in $output_dir."

  # Extract specific columns (1, 6, and 11) from the generated PAF files
  filtered_output1="${filtered_output_dir}/${sample_id}_hapdup_phased_1_filtered.paf"
  filtered_output2="${filtered_output_dir}/${sample_id}_hapdup_phased_2_filtered.paf"

  # Use awk to extract columns 1, 6, and 11 and save to new filtered output files
  echo "Extracting columns from $output1 to $filtered_output1..."
  awk '{print $1, $6, $11}' "$output1" > "$filtered_output1"

  echo "Extracting columns from $output2 to $filtered_output2..."
  awk '{print $1, $6, $11}' "$output2" > "$filtered_output2"

done

echo "Minimap2 processing and column extraction complete. Filtered outputs saved in $filtered_output_dir."

