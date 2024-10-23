#!/bin/bash
#PBS -l ncpus=8,mem=32GB,walltime=04:00:00,storage=gdata/te53+gdata/if89
#PBS -N bamextract
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/
#PBS -q normalbw

#Author: Sarah Jackson

set -e
set -o pipefail
set -u


#Pangenome VCF = /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/NEW/hprc-v1.1-mc-chm13-vcfbub.new.vcf.gz
#Jellyfish K-mer Counts = /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/NEW/kmercounts_chm13.jf
#Bed File of Loci = /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/NEW/blgrp_coords_chm13_anchorflank_named.bed
#List of Input files = /g/data/te53/sj2852/blgrp/analyses/reference_alignment_initial/NEW_alignments/combined_fastq

module load htslib/1.9
module load samtools
module load singularity

#SINGLE SAMPLES
#singularity exec /g/data/if89/singularityimg/locityper_0.15.2.sif /home/locityper/target/release/locityper genotype -i /g/data/te53/sj2852/blgrp/tmp/new_alignment/fastq_files_chm13/GM18501.blgrp.fastq.gz -d /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/database -p /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/bg -o /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/gts

#Following steps to be run embarrassingly parallel 

inputfastq=$INPUT_FQ
preproc_dir=$PREPROC_DIR
output_genotype=$OUTPUT_DIR
sid=$SAMPLEID


echo "inputfastq: $inputfastq"
echo "preproc_dir: $preproc_dir"
echo "output_genotype: $output_genotype"
echo "sid: $sid"

if [ -z "$OUTPUT_DIR" ]; then
    echo "ERROR: OUTPUT_DIR is not set"
    exit 1
fi

mkdir -p ${output_genotype}

# Check if the output directory has data
if [ -z "$(ls -A $output_genotype)" ]; then
    echo "Output directory is empty. Running genotyping."
    singularity exec /g/data/if89/singularityimg/locityper_0.15.2.sif /home/locityper/target/release/locityper genotype -i ${inputfastq} -d /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/NEW/database -p ${preproc_dir} -o ${output_genotype}
else
    echo "Genotyping output directory is not empty for sample $sid. Skipping genotyping."
fi



