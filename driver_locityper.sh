#!/bin/bash
#PBS -l ncpus=8,mem=32GB,walltime=04:00:00,storage=gdata/te53+gdata/if89
#PBS -N locityper_preproc_NCIG
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/
#PBS -q normalbw
#PBS -P te53

#Author: Sarah Jackson

set -e
set -o pipefail
set -u

module load htslib/1.9
module load samtools
module load singularity

#base_input_fastq_chm13="/g/data/te53/sj2852/blgrp/analyses/merge_alignment/kgp_chm13_fastq"
#base_input_fastq_grch38="/g/data/te53/sj2852/blgrp/analyses/merge_alignment/kgp_grch38_fastq"
base_output_preproc="/g/data/te53/sj2852/blgrp/analyses/locityper_kgp/NEW/bg"
base_output_genotype="/g/data/te53/sj2852/blgrp/analyses/locityper_kgp/NEW/gts"


#Run preprocessing
#for inputfastq in /g/data/te53/sj2852/blgrp/analyses/reference_alignment_initial/NEW_alignments/NCIG_all/ncig_fastq/*_pass.fastq.gz; 
    #do sid=$(basename "$inputfastq" _pass.fastq.gz); 
    #output_preproc="${base_output_preproc}/${sid}"

    #if [ ! -d "$output_preproc" ]; then

    #echo "starting preprocessing for sample: $sid, using input file: $inputfastq"

    #singularity exec /g/data/if89/singularityimg/locityper_0.15.2.sif /home/locityper/target/release/locityper preproc --tech nanopore -i "${inputfastq}" -r /g/data/te53/sj2852/blgrp/referenceresource/chm13-t2t.ebv.phix.chrQ.xy.fa -j /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/NEW/kmercounts_chm13.jf -o "${output_preproc}" 

    #qsub -P te53 -v INPUT_FQ="$inputfastq",OUTPUT_DIR="$output_preproc",SAMPLEID="$sid" /g/data/te53/sj2852/blgrp/scripts/alignment_assembly_code/locityper/locityper_preproc.sh; 
    #fi
#done

#Run genotyping
for inputfastq in /g/data/te53/sj2852/blgrp/analyses/reference_alignment_initial/NEW_alignments/combined_fastq/*_pass.fastq.gz;
    do
    sid=$(basename $inputfastq _pass.fastq.gz)
    output_preproc="${base_output_preproc}/${sid}"
    output_genotype="${base_output_genotype}/${sid}"

    if [ -d "${output_preproc}" ]; then

    mkdir -p ${output_genotype}

    qsub -P te53 -v INPUT_FQ=$inputfastq,PREPROC_DIR=$output_preproc,OUTPUT_DIR=$output_genotype,SAMPLEID=$sid /g/data/te53/sj2852/blgrp/scripts/alignment_assembly_code/locityper/locityper_genotype.sh
    fi
done