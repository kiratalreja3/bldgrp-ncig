#!/bin/bash
#PBS -q copyq
#PBS -P te53
#PBS -l ncpus=1,walltime=10:00:00,mem=4GB,storage=gdata/te53+gdata/if89
#PBS -N dl1kgp
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/

set -ex
module use -a /g/data/if89/apps/modulefiles
module load samtools aws-cli


s3baseurl="1000g-ont/ALIGNMENT_AND_ASSEMBLY_DATA/FIRST_100/IN-HOUSE_MINIMAP2/CHM13/"
blgrpchm13bed="/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/NEW_COORDS/chm13/blgrp_coords_chm13_anchorflank.bed"
#"/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/blgrp.coords.chm13.mergeflank.sorted.fix.chr17added.bed"
workdir=/g/data/te53/sj2852/blgrp/data/NEW_KGP_DATA/onekgpchm13_bam

cd ${workdir}
dllogfile=onekgp_ont_chm13bam_blgrp_download.log

for i in `aws s3 ls "s3://${s3baseurl}" --no-sign-request`
do
	if [ $i != "PRE" ]
	then
		sid=$(echo $i | cut -f1 -d '-')
		aws s3 ls "s3://${s3baseurl}${i}" --no-sign-request >$sid.filelist
		fname=$(grep -P ".bam$" ${sid}.filelist | sed -r 's/.*\s+//')

		#If output already exists, skip
		 if [ -f "${sid}.blgrp.bam" ]; then
            echo "Output file ${sid}.blgrp.bam already exists. Skipping..."
            continue  
        fi

		samtools view -h --region-file ${blgrpchm13bed} --bam --output ${sid}.blgrp.bam##idx##${sid}.blgrp.bam.bai --write-index  https://s3.amazonaws.com/${s3baseurl}${i}${fname}
		echo ${sid} ${fname} $(date) >>${dllogfile}
		rm ${sid}.filelist ${fname}.bai
	fi
done

