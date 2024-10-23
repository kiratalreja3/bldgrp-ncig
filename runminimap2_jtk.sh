#!/bin/bash
#PBS -l ncpus=8,mem=32GB,walltime=06:00:00,storage=gdata/te53+gdata/if89
#PBS -N minimap2_hapdup
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/
#PBS -q normal
#PBS -P te53

#File Locations
#READS=
#REFERENCE="/g/data/te53/sj2852/blgrp/referenceresource/chm13-t2t.ebv.phix.chrQ.xy.fa"
#REGIONS="/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/chm13/blgrp.coords.chm13.mergeflank.sorted.fix.chr17added.bed"
#CONFIG=""
#OUTPUT_DIR=""

set -e
set -o pipefail
set -u

usage() {
        echo "Usage: qsub -v READS=input.fasta,REFERENCE=chm13-t2t.ebv.phix.chrQ.xy.fa,REGIONS=blgrp.coords.chm13.mergeflank.sorted.fix.chr17added.bed,CONFIG= OUTPUT_DIR=outdir,SAMPLEID=sampleid /g/data/te53/sj2852/blgrp/scripts/run_and_submit_scripts/runminimap2_hapdup.sh" >&2
        echo
        exit 1
}

#Check for Required Files
if [ -z "${READS}" ] || [ -z "${REFERENCE}" ] || [ -z "${REGION}" ]|| [ -z "${CONFIG}" ]; then
    echo "Missing required variables" >&2
    exit 1
fi

#Check for Output File
output="${OUTPUT_DIR}/${SAMPLEID}.out.gfa"

if [ -f "$output" ]; then
    echo "Output file ${output} already exists. Skipping this sample."
    exit 0
fi

module load minimap2/2.28
module load samtools

echo "Generating Minimap Alignment"


echo "Extracting Reads for Target Ranges"

#edit
minimap2 -x map-ont -t $THREADS --secondary=no -a $REFERENCE $READS |
    samtools sort -@$THREADS -OBAM > aln.bam
samtools index aln.bam
samtools view -OBAM aln.bam $REGION |
    samtools fasta > reads.fasta