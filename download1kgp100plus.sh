#!/bin/bash
#PBS -q copyq
#PBS -P te53
#PBS -l ncpus=1,walltime=10:00:00,mem=4GB,storage=gdata/te53+gdata/if89
#PBS -N dl1kgp_100plus
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/

set -ex
module use -a /g/data/if89/apps/modulefiles
module load samtools aws-cli

#https://s3.amazonaws.com/1000g-ont/ALIGNMENT_AND_ASSEMBLY_DATA/100_PLUS/IN-HOUSE_MINIMAP2/GM18507-ONT-hg38-R9-LSK110-dorado050_sup_5mCG_v33/GM18507-ONT-hg38-R9-LSK110-dorado050_sup_5mCG_v33.phased.bam
s3baseurl="1000g-ont/ALIGNMENT_AND_ASSEMBLY_DATA/100_PLUS/IN-HOUSE_MINIMAP2/"
blgrpgrch38bed="/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/NEW_COORDS/grch38/blgrp_coords_grch38_anchorflank.bed"
workdir=/g/data/te53/sj2852/blgrp/data/NEW_KGP_DATA/onekgp100plusgrch38_bam

cd ${workdir}
dllogfile=onekgp_ont_grch38bam_blgrp_download.log

for i in `aws s3 ls "s3://${s3baseurl}" --no-sign-request`
do
	if [ $i != "PRE" ]
	then
		sid=$(echo $i | cut -f1 -d '-')
		aws s3 ls "s3://${s3baseurl}${i}" --no-sign-request >$sid.filelist
		fname=$(grep -P ".bam$" ${sid}.filelist | sed -r 's/.*\s+//')

		#Define output file paths
		        output_bam="${sid}.blgrp.bam"
        		output_bai="${sid}.blgrp.bam.bai"

		#if file already exists, skip it
		        if [ -f "$output_bam" ] && [ -f "$output_bai" ]; then
			            echo "Skipping ${sid}, output files already exist."
						continue
				fi
	
		#processing file
		echo "Processing ${sid}..."
		samtools view -h --region-file ${blgrpgrch38bed} --bam --output ${sid}.blgrp.bam##idx##${sid}.blgrp.bam.bai --write-index  https://s3.amazonaws.com/${s3baseurl}${i}${fname}
		echo ${sid} ${fname} $(date) >>${dllogfile}
		rm ${sid}.filelist
	fi
done

