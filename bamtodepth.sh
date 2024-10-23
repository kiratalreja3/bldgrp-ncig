#!/bin/bash
#PBS -N bedcov
#PBS -P te53
#PBS -q normal
#PBS -l ncpus=48,mem=192GB,walltime=01:00:00,storage=gdata/te53+gdata/if89
 
#this file calculates the read depth from my KGP-ONT files for each blood group gene, the depth is normalised for each blood group gene
#submit using the submitbamtodepth.sh script (/g/data/te53/sj2852/blgrp/scripts/run_and_submit_scripts/submitbamtodepth.sh)
 
module load samtools
module load parallel


bambase="$(basename ${bam} .bam)"
 
# Function to calculate depth for each BED region
calculate_depth() {
    region=$1
    samtools bedcov <(echo "$region") ${bam} | awk '{window=$3-$2; print $1, $2, $3, $4/window}'
}
 
export -f calculate_depth
 
# Run depth calculation in parallel
parallel --will-cite -a ${bed} -j ${PBS_NCPUS} calculate_depth > "${outdir}/${bambase}.depth.tmp.bed"
sort -k1,1 -k2,2n "${outdir}/${bambase}.depth.tmp.bed" > "${outdir}/${bambase}.depth.bed"
rm "${outdir}/${bambase}.depth.tmp.bed"