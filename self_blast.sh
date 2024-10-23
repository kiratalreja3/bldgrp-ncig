#!/bin/bash
#PBS -l ncpus=8,mem=32GB,walltime=03:00:00,storage=gdata/te53+gdata/if89
#PBS -N blast_table
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/
#PBS -q normal
#PBS -P te53

#Author: Hardip Patel

#this script takes the MB fasta file of a gene, creates a database and then aligns to itself to identify regions of overlap

module load blast/2.11.0

#First Create the Database
makeblastdb -dbtype nucl -in /g/data/te53/sj2852/blgrp/analyses/anchors/GYPA_GYPB_anchor_identifying/300kb_5_100kb_3/gypa_gypb_300_100kb.fa

#Then make the table
blastn -query /g/data/te53/sj2852/blgrp/analyses/anchors/GYPA_GYPB_anchor_identifying/300kb_5_100kb_3/gypa_gypb_300_100kb.fa -db /g/data/te53/sj2852/blgrp/analyses/anchors/GYPA_GYPB_anchor_identifying/300kb_5_100kb_3/gypa_gypb_300_100kb.fa -outfmt 6 -out /g/data/te53/sj2852/blgrp/analyses/anchors/GYPA_GYPB_anchor_identifying/300kb_5_100kb_3/gypa.vs.gypb.300_100kb.tbl