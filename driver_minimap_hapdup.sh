#!/bin/bash

for file in /g/data/te53/sj2852/blgrp/analyses/flye_assemblies/NEW_FLYE/hapdup/assembly_fasta/*.flyeasm.fasta; 
    do sid=$(basename $file .flyeasm.fasta); 
    input_fastq1=/g/data/te53/sj2852/blgrp/analyses/reference_alignment_initial/NEW_alignments/NCIG_all/ncig_fastq/${sid}_pass.fastq.gz
    #input_fastq2=/g/data/te53/sj2852/blgrp/analyses/flye_assemblies/NEW_FLYE/hapdup/fastq_files/${sid}.blgrp.fastq.gz
    outdir=/g/data/te53/sj2852/blgrp/analyses/flye_assemblies/NEW_FLYE/hapdup/lr_mapping_bam

    if [[ -f "$input_fastq1" ]]; then
        input_fastq="$input_fastq1"
    #elif [[ -f "$input_fastq2" ]]; then
        #input_fastq="$input_fastq2"
    else
        echo "No input FASTQ files found for $sid. Skipping..."
        continue
    fi

    qsub -v ASSEMBLY_FASTA=$file,READS_FASTQ=$input_fastq,OUTPUT_DIR=$outdir,SAMPLEID=$sid /g/data/te53/sj2852/blgrp/scripts/run_and_submit_scripts/runminimap2_hapdup.sh; 
done 
