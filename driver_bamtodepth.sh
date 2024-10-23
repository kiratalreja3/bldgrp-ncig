#!/bin/bash

bed=/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/NEW_COORDS/chm13/blgrp_coords_chm13_anchorflank.bed
outdir=/g/data/te53/sj2852/blgrp/metadata/assembly_stats/NEW/bedcov_NCIG


for file in /g/data/te53/sj2852/blgrp/analyses/reference_alignment_initial/NEW_alignments/NCIG_all/ncig_bam_and_bai/*.bam; do
    # Check if the current file is .bam NOT .bam.bai
     if [[ -f "$file" && ! "$file" =~ \.bam\.bai$ ]]; then
        qsub -v bam=$file,bed=$bed,outdir=$outdir /g/data/te53/sj2852/blgrp/scripts/stats_by_blgrp/bamtodepth.sh
    fi

 done   
