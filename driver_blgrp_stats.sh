#!/bin/bash

for file in /g/data/te53/sj2852/blgrp/data/KGPdata/onekgpchm13_bam/*.bam;
    do sid=$(basename $file .bam);
    blgrp_regions=/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/bloodcoords_withregion_chm13.bed;
    outdir=/g/data/te53/sj2852/blgrp/metadata/assembly_stats/by_blgrp_region/chm13_region_stats;

    qsub -P te53 -v INPUT_BAM=$file,INPUT_REF=$blgrp_regions,OUTPUT_DIR=$outdir,SAMPLEID=$sid /g/data/te53/sj2852/blgrp/scripts/blgrp_stats.sh;
done



