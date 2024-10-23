#!/bin/bash


for file in /g/data/te53/sj2852/blgrp/analyses/flye_assemblies/NEW_FLYE/*.flyeasm.fasta; 
    do sid=$(basename $file .blgrp.fasta); 
    ref_file=/g/data/te53/sj2852/blgrp/referenceresource/chm13-t2t.ebv.phix.chrQ.xy.fa
    outdir=/g/data/te53/sj2852/blgrp/analyses/flye_assemblies/NEW_FLYE/fasta_output

    qsub -v INPUT_FASTA=$file,INPUT_REF=$ref_file,OUTPUT_DIR=$outdir,SAMPLEID=$sid /g/data/te53/sj2852/blgrp/scripts/run_and_submit_scripts/runminimap2_flye.sh; 
done 
