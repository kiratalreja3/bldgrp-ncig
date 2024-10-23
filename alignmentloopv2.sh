#!/bin/bash
#PBS -l ncpus=8,mem=32GB,walltime=06:00:00,storage=gdata/te53+gdata/if89
#PBS -N bamextract_NCIG
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/
#PBS -q normalbw

#SINGLE SAMPLES
#bed_file="/g/data/te53/sj2852/blgrp/metadata/blgrp.coords.chm13.mergeflank.sorted.bed"
#output_folder="/g/data/te53/sj2852/blgrp/analyses/merge_alignment/"
#inputbam="/g/data/te53/sj2852/bam_files/PGXX22494_pass.bam"

#module load samtools
#echo "Bed file: $bed_file"
#while IFS=$'\t' read -r chromosome start end; do
    #output_file="$output_folder/${chromosome}_${start}-${end}_aligned_reads.bam"
   
#output_file="$output_folder/$(basename ${inputbam} .bam).blgrpsubset.bam"

#echo "Generating BAM"
#samtools view --threads ${PBS_NCPUS} -h -M -L ${bed_file} --write-index ${inputbam} -b -o ${output_file}##idx##${output_file}.bai

#echo "Generating Fastq"
#samtools fastq ${output_file} | gzip > ${output_folder}/$(basename ${inputbam} .bam).blgrpsubset.fastq.gz
 
#flye
#    echo "Aligned reads for $chromosome:$start-$end saved to $output_file"


#done < "$bed_file"


#MULTIPLE SAMPLES

bed_file="/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/NEW_COORDS/chm13/blgrp_coords_chm13_anchorflank_named.bed" #for chm13
#bed_file="/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/blgrp.coords.grch38.mergeflank.sorted.chr17added.bed" #for grch38

#output_folder="/g/data/te53/sj2852/blgrp/analyses/merge_alignment/1kgp"
output_folder="/g/data/te53/sj2852/blgrp/analyses/reference_alignment_initial/NEW_alignments/kgp_chm13_all/kgp_chm13_bam_and_bai"

#input_folder="/g/data/te53/sj2852/blgrp/data/NEW_KGP_DATA/onekgpchm13_bam" #for chm13
#input_folder="/g/data/te53/sj2852/blgrp/data/KGPdata/onekgp100plusgrch38_bam"  #for grch38
input_folder="/g/data/te53/sj2852/blgrp/data/NEW_KGP_DATA/onekgpchm13_bam"

module load samtools

for inputbam in "${input_folder}"/*.bam; do
     #Check if the current item in the loop is a regular file
    if [ -f "$inputbam" ]; then
        #Define output file paths
        output_file="${output_folder}/$(basename ${inputbam} .bam).blgrpsubset.bam"
        output_fastq="${output_folder}/$(basename ${inputbam} .bam).fastq.gz"
        
        #echo "Generating BAM: ${inputbam}"
        #samtools view --threads "${PBS_NCPUS}" -h -M -L "${bed_file}" --write-index "${inputbam}" -b -o "${output_file}##idx##${output_file}.bai"
        
        if [ -f "$output_file" ]; then
            echo "Generating Fastq: ${output_fastq}"
        
        samtools fastq "${output_file}" | gzip > "${output_fastq}"

        else   
            echo "Warning: BAM file '$output_file' does not exist. Skipping."
        fi
    else   
        echo "Warning: Input BAM file '$inputbam' does not exist."
    fi
done