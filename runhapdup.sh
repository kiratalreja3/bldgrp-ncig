#!/bin/bash
#PBS -l ncpus=8,mem=32GB,walltime=3:00:00,storage=gdata/te53+gdata/if89
#PBS -N runhapdup
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/
#PBS -q normal
#PBS -P te53

#File Locations
#INPUT_FASTA="/g/data/te53/sj2852/blgrp/analyses/flye_assemblies/NEW_FLYE/hapdup/assembly_fasta"
#INPUT_BAM="/g/data/te53/sj2852/blgrp/analyses/flye_assemblies/NEW_FLYE/hapdup/lr_mapping_bam*.alignment.bam"
#OUTPUT_DIR="/g/data/te53/sj2852/blgrp/analyses/flye_assemblies/NEW_FLYE/hapdup"


inputfasta=$INPUT_FASTA
inputbam=$INPUT_BAM
outdir=$OUTPUT_DIR
sid=$SAMPLEID


echo "inputfasta: $inputfasta"
echo "inputbam: $inputbam"
echo "outdir: $outdir"
echo "sid: $sid"

set -ex
set -o pipefail
set -u

usage() {
        echo "Usage: qsub -v INPUT_FASTA=input.flyeasm.fasta,INPUT_BAM=input.alignment.bam,OUTPUT_DIR=outdir,SAMPLEID=sampleid" /g/data/te53/sj2852/blgrp/scripts/run_and_submit_scripts/submithapdup.sh >&2
        echo
        exit 1
}

#Check for Required Files
if [ -z "${INPUT_FASTA}" ] || [ -z "${INPUT_BAM}" ] || [ -z "${OUTPUT_DIR}" ]; then
    echo "Missing required variables" >&2
    exit 1
fi

module load singularity
module load samtools

echo "Generating HapDup Assembly"

singularity exec \
    --bind /g/data/te53/sj2852/blgrp/analyses/flye_assemblies/NEW_FLYE/hapdup/lr_mapping_bam \
    --bind /g/data/te53/sj2852/blgrp/analyses/flye_assemblies/NEW_FLYE/hapdup/assembly_fasta \
    /g/data/if89/singularityimg/hapdup_0.12.sif hapdup \
    --assembly ${inputfasta} \
    --bam ${inputbam} \
    --out-dir ${outdir} \
    -t ${PBS_NCPUS} \
    --rtype ont

#Single File Test

#HD_DIR="/g/data/te53/sj2852/blgrp/tmp/hapdup_test"

#singularity exec /g/data/if89/singularityimg/hapdup_0.12.sif hapdup \
    #--assembly $HD_DIR/GM18501.flyeasm.fasta \
    #--bam $HD_DIR/GM18501.alignment.bam \
    #--out-dir $HD_DIR/hapdup \
    #-t ${PBS_NCPUS} \
    #--rtype ont 