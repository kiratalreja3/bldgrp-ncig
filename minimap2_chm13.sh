#!/bin/bash
#PBS -l ncpus=48,mem=192GB,walltime=01:00:00,storage=gdata/te53+gdata/if89
#PBS -N bamextract
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/
#PBS -q normal

#inputreference="/g/data/te53/humanreference/chm13-t2t/minimapidx/chm13-t2t.ebv.phix.chrQ.xy.fa.mm2idx"
#inputassembly="/g/data/te53/sj2852/blgrp/analyses/flye_assemblies/fasta_output_chm13/GM18501.flyeasm.fasta"
#outdir="/g/data/te53/sj2852/blgrp/analyses/minimap/1kgp_minimap"

module load minimap2/2.28


#Run Minimap2 (MULTIPLE)
minimap2 -x asm20 --cs --secondary=no -t ${PBS_NCPUS} "${inputreference}" "${inputassembly}" > ${outdir}/alignment.paf

#Single File
#minimap2 -ax map-ont /g/data/te53/sj2852/blgrp/referenceresource/chm13-t2t.ebv.phix.chrQ.xy.fa /g/data/te53/sj2852/blgrp/tmp/new_alignment/chm13_fastq/GM18501.blgrp.fastq.gz > /g/data/te53/sj2852/blgrp/tmp/new_alignment/GM18501.aln.sam

#minimap2 --cs /g/data/te53/sj2852/blgrp/referenceresource/chm13-t2t.ebv.phix.chrQ.xy.fa /g/data/te53/sj2852/blgrp/analyses/anchors/test_duplicate_coords.fa > /g/data/te53/sj2852/blgrp/analyses/200bp_investigation/200bp_alignment_global_assemblies/chm13_new.aln.paf

#minimap2 --cs  /g/data/te53/sj2852/blgrp/referenceresource/GRCh38_full_plus_hs38d1_IMGTHLAv3.54.0_IPDKITv2.13.0_patchesv14_ctrlseq_masked.xy.fa /g/data/te53/sj2852/blgrp/metadata/blgrp_coords/grch38/bloodcoords_200bpflanks_grch38.fa > /g/data/te53/sj2852/blgrp/analyses/200bp_investigation/200bp_alignment_global_assemblies/grch38.aln.paf

#minimap2 --cs /g/data/te53/sj2852/blgrp/referenceresource/chm13-t2t.ebv.phix.chrQ.xy.fa /g/data/te53/sj2852/blgrp/analyses/anchors/full_script_test/blgrp_coords.fasta > /g/data/te53/sj2852/blgrp/analyses/anchors/full_script_test/anchor_alignment_initial.paf