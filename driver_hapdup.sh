#!/bin/bash

for file in /g/data/te53/sj2852/blgrp/analyses/flye_assemblies/NEW_FLYE/hapdup/assembly_fasta/*.flyeasm.fasta; 
    do sid=$(basename $file .flyeasm.fasta); 

    #define paths
    inputbam=/g/data/te53/sj2852/blgrp/analyses/flye_assemblies/NEW_FLYE/hapdup/lr_mapping_bam/${sid}.alignment.bam
    outdir=/g/data/te53/sj2852/blgrp/analyses/flye_assemblies/hapdup/${sid}

    mkdir -p $outdir

    qsub -v INPUT_FASTA=$file,INPUT_BAM=$inputbam,OUTPUT_DIR=$outdir,SAMPLEID=$sid /g/data/te53/sj2852/blgrp/scripts/run_and_submit_scripts/runhapdup.sh;
    
done 
