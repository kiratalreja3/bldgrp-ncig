#!/bin/bash


for file in /g/data/te53/sj2852/blgrp/analyses/flye_assemblies/fasta_output_all/*.flyeasm.fasta; 
    do sid=$(basename $file .flyeasm.fasta); 
    #input_fastq=/g/data/te53/sj2852/blgrp/analyses/flye_assemblies/fastq_files_for_hapdup/${sid}.blgrp.fastq
    ref_file=/g/data/te53/sj2852/blgrp/referenceresource/chm13-t2t.ebv.phix.chrQ.xy.fa
    outdir=/g/data/te53/sj2852/blgrp/analyses/minimap_flye_assemblies/flye_align_to_chm13

    qsub -v INPUT_FASTA=$file,INPUT_REF=$ref_file,OUTPUT_DIR=$outdir,SAMPLEID=$sid /g/data/te53/sj2852/blgrp/scripts/run_and_submit_scripts/runminimap2.sh; 
done 
