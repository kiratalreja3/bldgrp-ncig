#!/bin/bash
#PBS -l ncpus=12,mem=120GB,walltime=04:00:00,storage=gdata/te53+gdata/if89
#PBS -N jellyfish
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/
#PBS -q normalbw

#Original Authors : Hardip Patel, Kirat Alreja, Arthur Georges
#Modified by: Sarah Jackson

input_file="/g/data/te53/sj2852/blgrp/referenceresource/chm13-t2t.ebv.phix.chrQ.xy.fa"
output_dir="/g/data/te53/sj2852/blgrp/analyses/locityper_kgp/NEW"

set -ex

module load jellyfish/2.3.0 utils/0.0


jellyfish count --canonical --lower-count 2 --out-counter-len 2 --mer-len 25 \
    --threads=${PBS_NCPUS} --size 3G --output "${output_dir}/kmercounts_chm13.jf" "${input_file}"

