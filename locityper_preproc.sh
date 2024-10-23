#!/bin/bash
#PBS -l ncpus=4,mem=16GB,walltime=04:00:00,storage=gdata/te53+gdata/if89
#PBS -N locityper_preproc
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/
#PBS -q normalbw

#Author: Sarah Jackson

set -e
set -o pipefail
set -u

#Pangenome Raw = /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/NEW/hprc-v1.1-mc-chm13.raw.new.vcf.gz
#Pangenome VCF = /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/NEW/hprc-v1.1-mc-chm13-vcfbub.new.vcf.gz
#Jellyfish K-mer Counts = /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/kmercounts_chm13_new.jf
#Bed File of Loci = /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/bloodcoords_withregion.bed
#List of Input files = /g/data/te53/sj2852/blgrp/analyses/merge_alignment/combined_fastq

module load htslib/1.9
module load samtools
module load singularity

#Transform pangenome raw VCF so that no variants overlap each other
#singularity exec /g/data/if89/singularityimg/pggb_latest.sif vcfbub -l 0 -i /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/NEW/hprc-v1.1-mc-chm13.raw.new.vcf.gz | bgzip > /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/hprc-v1.1-mc-chm13-vcfbub.new1.vcf.gz

#index the pangenome VCF
#tabix -p vcf /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/NEW/hprc-v1.1-mc-chm13-vcfbub.new.vcf.gz


# locityper add command 
#singularity exec /g/data/if89/singularityimg/locityper_0.15.2.sif /home/locityper/target/release/locityper add -e 3000000 -d /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/NEW/database -v /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/NEW/hprc-v1.1-mc-chm13-vcfbub.new.vcf.gz -r /g/data/te53/sj2852/blgrp/referenceresource/chm13-t2t.ebv.phix.chrQ.xy.fa -j /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/NEW/kmercounts_chm13.jf -L /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/NEW/blgrp_coords_chm13_anchorflank_named.bed


#run preprocessing
singularity exec /g/data/if89/singularityimg/locityper_0.15.2.sif /home/locityper/target/release/locityper preproc --tech nanopore -i "${inputfastq}" -r /g/data/te53/sj2852/blgrp/referenceresource/chm13-t2t.ebv.phix.chrQ.xy.fa -j /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/NEW/kmercounts_chm13.jf -o "${output_preproc}" 




#SINGLE SAMPLES

#run WGS dataset preprocessing TEST
#singularity exec /g/data/if89/singularityimg/locityper_0.15.2.sif /home/locityper/target/release/locityper preproc --tech nanopore -i /g/data/te53/sj2852/blgrp/tmp/new_alignment/fastq_files_chm13/GM18501.blgrp.fastq.gz -r /g/data/te53/sj2852/blgrp/referenceresource/chm13-t2t.ebv.phix.chrQ.xy.fa -j /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/kmercounts_chm13_new.jf -o /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/bg -b chr7:30984043-31067163

#targeted genotyping
#singularity exec /g/data/if89/singularityimg/locityper_0.15.2.sif /home/locityper/target/release/locityper genotype -i /g/data/te53/sj2852/blgrp/tmp/new_alignment/fastq_files_chm13/HG01395.blgrp.fastq.gz -d /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/database -p /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/bg -o /g/data/te53/sj2852/blgrp/analyses/locityper_kgp/gts

