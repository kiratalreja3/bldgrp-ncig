#!/bin/bash


for file in /g/data/te53/sj2852/blgrp/data/KGPdata/onekgpchm13_bam/*.blgrp.bam; 
    do sid=$(basename $file .blgrp.bam); 
    outdir=/g/data/te53/sj2852/blgrp/analyses/flye_assemblies/fastq_files_for_hapdup

    qsub -P te53 -v input_bam=$file,output_folder=$outdir,SAMPLEID=$sid /g/data/te53/sj2852/blgrp/scripts/alignment_assembly_code/bamtofastq_copy.sh; 
done 
