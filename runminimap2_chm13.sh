#!/bin/bash
#PBS -l ncpus=8,mem=32GB,walltime=10:00:00,storage=gdata/te53+gdata/if89
#PBS -N minimap_chm13
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/
#PBS -q normal
#PBS -P te53

#File Locations
#INPUT_FASTA= #chm13.fasta files 
#INPUT_REF="/g/data/te53/sj2852/blgrp/referenceresource/chm13-t2t.ebv.phix.chrQ.xy.fa"
#OUTPUT_DIR="/g/data/te53/sj2852/blgrp/analyses/reference_alignment_initial/NEW_alignments/kgp_chm13"


set -e
set -o pipefail
set -u

usage() {
        echo "Usage: qsub -v INPUT_FASTA=input.flyeasm.fasta,INPUT_REF=input.ebv.phix.chrQ.xy.fa.mm2idx,OUTPUT_DIR=outdir,SAMPLEID=sampleid /g/data/te53/sj2852/blgrp/scripts/run_and_submit_scripts/runminimap2_chm13.sh" >&2
        echo
        exit 1
}

#Check for Required Files
if [ -z "${INPUT_FASTA}" ] || [ -z "${INPUT_REF}" ] || [ -z "${OUTPUT_DIR}" ]; then
    echo "Missing required variables" >&2
    exit 1
fi

#Check for Output File
output_bam="${OUTPUT_DIR}/${SAMPLEID}.alignment.bam"

if [ -f "$output_bam" ]; then
    echo "Output file ${output_bam} already exists. Skipping this sample."
    exit 0
fi

module load minimap2/2.28
module load samtools

echo "Generating Minimap Alignment"
minimap2 -x asm20 --cs --secondary=no -t ${PBS_NCPUS} "${INPUT_REF}" "${INPUT_FASTA}" > ${OUTPUT_DIR}/${SAMPLEID}.alignment.paf


#.alignment.paf or .alignment.sam depending on what you want to do with the output...