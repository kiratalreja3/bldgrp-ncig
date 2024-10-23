#!/bin/bash
#PBS -l ncpus=48,mem=192GB,walltime=02:00:00,storage=gdata/te53+gdata/if89
#PBS -N minimap2_hapdup
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/
#PBS -q normal
#PBS -P te53

#File Locations
#ASSEMBLY_FASTA="/g/data/te53/sj2852/blgrp/analyses/flye_assemblies/NEW_FLYE/hapdup/assembly_fasta/*.flyeasm.fasta"
#READS_FASTQ="/g/data/te53/sj2852/blgrp/analyses/flye_assemblies/NEW_FLYE/hapdup/fastq_files/*.blgrp.fastq"
#OUTPUT_DIR="/g/data/te53/sj2852/blgrp/analyses/flye_assemblies/NEW_FLYE/hapdup/lr_mapping_bam"


set -e
set -o pipefail
set -u

usage() {
        echo "Usage: qsub -v ASSEMBLY_FASTA=input.flyeasm.fasta,READS_FASTQ=input.blgrp.fastq,OUTPUT_DIR=outdir,SAMPLEID=sampleid /g/data/te53/sj2852/blgrp/scripts/run_and_submit_scripts/runminimap2_hapdup.sh" >&2
        echo
        exit 1
}

#Check for Required Files
if [ -z "${ASSEMBLY_FASTA}" ] || [ -z "${READS_FASTQ}" ] || [ -z "${OUTPUT_DIR}" ]; then
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
minimap2 -ax map-ont -t ${PBS_NCPUS} "${ASSEMBLY_FASTA}" "${READS_FASTQ}" | samtools sort -@ ${PBS_NCPUS} -m 4G > "${OUTPUT_DIR}"/"${SAMPLEID}".alignment.bam

echo "Indexing BAM file"
samtools index -@ ${PBS_NCPUS} "${OUTPUT_DIR}"/"${SAMPLEID}".alignment.bam

#Single Files

#echo "Generating Minimap Alignment"
#minimap2 -ax map-ont -t ${PBS_NCPUS} /g/data/te53/sj2852/blgrp/analyses/flye_assemblies/NEW_FLYE/NCIG_flye/PGXX22492.flyeasm/PGXX22492.flyeasm.fasta /g/data/te53/sj2852/blgrp/analyses/flye_assemblies/NEW_FLYE/hapdup/fastq_files/PGXX22492_pass.fastq.gz | samtools sort -@ ${PBS_NCPUS} -m 4G > /g/data/te53/sj2852/blgrp/analyses/flye_assemblies/NEW_FLYE/hapdup/lr_mapping_bam/PGXX22492.alignment.bam

#echo "Indexing BAM file"
#samtools index -@ ${PBS_NCPUS} /g/data/te53/sj2852/blgrp/analyses/flye_assemblies/NEW_FLYE/hapdup/lr_mapping_bam/PGXX22492.alignment.bam
