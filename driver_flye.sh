#!/bin/bash


for file in /g/data/te53/sj2852/blgrp/analyses/reference_alignment_initial/NEW_alignments/kgp_chm13_all/kgp_chm13_fastq/*.blgrp.fastq.gz; 
    do sid=$(basename $file .blgrp.fastq.gz); 
    outdir=/g/data/te53/sj2852/blgrp/analyses/flye_assemblies/NEW_FLYE/chm13_flye/$sid.flyeasm; 

    qsub -v INPUT_FQ=$file,OUTPUT_DIR=$outdir,SAMPLEID=$sid /g/data/te53/sj2852/blgrp/scripts/run_and_submit_scripts/flye/runflye.sh; 
done 


#modify file location and submit