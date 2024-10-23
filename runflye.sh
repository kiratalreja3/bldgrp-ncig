#!/bin/bash
#PBS -l ncpus=4,mem=32GB,walltime=03:00:00,storage=gdata/te53+gdata/if89
#PBS -N flye_grch38
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/
#PBS -q normal
#PBS -P te53

set -e
set -o pipefail
set -u
 
usage() {
        echo "Usage: qsub -v INPUT_FQ=input.fq.gz,OUTPUT_DIR=outputdir,SAMPLEID=sampleid /g/data/te53/sj2852/blgrp/scripts/alignment_assembly_code/runflye.sh'" >&2
        echo
        exit 1
}

module load Flye/2.9.3

if [[ -e ${OUTPUT_DIR}/${SAMPLEID}.flyeasm.done ]]
then
    echo "Nothing to be done"
    exit 0
fi

mkdir ${OUTPUT_DIR}

echo STARTTIME = $(date) > ${OUTPUT_DIR}/${SAMPLEID}.flyeasm.running
#echo "Generating Flye Assembly"
flye --nano-hq ${INPUT_FQ} --out-dir ${OUTPUT_DIR} --iterations 3 --threads ${PBS_NCPUS} 

rm -r ${OUTPUT_DIR}/00-assembly
rm -r ${OUTPUT_DIR}/10-consensus
rm -r ${OUTPUT_DIR}/20-repeat
rm -r ${OUTPUT_DIR}/30-contigger
rm -r ${OUTPUT_DIR}/40-polishing
rm ${OUTPUT_DIR}/assembly_graph.gv

sed "s/>/>${SAMPLEID}_/" ${OUTPUT_DIR}/assembly.fasta > $OUTPUT_DIR/${SAMPLEID}.flyeasm.fasta
mv ${OUTPUT_DIR}/assembly_graph.gfa $OUTPUT_DIR/${SAMPLEID}.flyeasm.gfa
mv ${OUTPUT_DIR}/assembly_info.txt $OUTPUT_DIR/${SAMPLEID}.flyeasm.info.txt

rm ${OUTPUT_DIR}/assembly.fasta

mv ${OUTPUT_DIR}/${SAMPLEID}.flyeasm.running ${OUTPUT_DIR}/${SAMPLEID}.flyeasm.done
echo ENDTIME = `date` >> ${OUTPUT_DIR}/${SAMPLEID}.flyeasm.done

#first run without the 2000 min overlap as a baseline, add this in afterwards to see what happens

#for file in /g/data/te53/sj2852/blgrp/analyses/merge_alignment/1kgp_fastq/*.fastq.gz; do sid=$(basename $file .blgrp.fastq.gz); outdir=/g/data/te53/sj2852/blgrp/analyses/flye_assemblies/$sid.flyeasm; echo qsub -v INPUT_FQ=$file,OUTPUT_DIR=$outdir,SAMPLEID=$sid  /g/data/te53/sj2852/blgrp/scripts/alignment_assembly_code/runflye.sh; done

