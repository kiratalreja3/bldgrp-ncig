#!/bin/bash
#PBS -l ncpus=8,mem=32GB,walltime=03:00:00,storage=gdata/te53+gdata/if89
#PBS -N bamextract
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/
#PBS -q normal
#PBS -P te53


module load samtools 

for bam in /g/data/te53/sj2852/blgrp/analyses/minimap_flye_assemblies/minimap_bam_new/*.bam; do     

samtools index $bam; done