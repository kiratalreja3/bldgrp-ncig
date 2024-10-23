#!/bin/bash
#PBS -l ncpus=4,mem=16GB,walltime=2:00:00,storage=gdata/te53+gdata/if89
#PBS -N full_anchor
#PBS -j oe
#PBS -o /g/data/te53/sj2852/blgrp/logs/
#PBS -q normal
#PBS -P te53

#Author: Sarah Jackson

#This script is the full pipeline from input coordinates to identifying anchor regions 

module load Rlib/4.3.1
module load samtools
module load minimap2/2.26
module load blast/2.14.1


#Input Files
input_bed_file="/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/chm13/bloodcoords_liftover_pre_mergeflank.bed" #Input bed should be in structure (chr start end gene) all tab separated
reference_fasta="/g/data/te53/sj2852/blgrp/referenceresource/chm13-t2t.ebv.phix.chrQ.xy.fa"
output_directory="/g/data/te53/sj2852/blgrp/tmp/full_coordinate_test/take_3"

#Output Directories
#mkdir -p "$output_directory"
#mkdir -p "$blast_output_dir"

#Different Output Paths
coordinate_anchors_bed="${output_directory}/coordinate_anchors.bed"
coordinate_anchors_cleaned_bed="${output_directory}/coordinate_anchors_cleaned.bed"
full_region_bed="${output_directory}/full_region.bed"
blast_reference_db="${output_directory}/blast_reference_db"
minimap_anchor_alignment_initial_paf="${output_directory}/minimap_anchor_alignment_initial.paf"
blast_anchor_alignment_initial_tbl="${output_directory}/blast_anchor_alignment_initial.tbl"
anchor_top3_tbl="${output_directory}/anchor_top_3.tbl"
anchor_filtered_tbl="${output_directory}/anchor_filtered.tbl"
blgrp_anchors_fasta="${output_directory}/blgrp_anchors.fasta"
labelled_minimap_alignment_paf="${output_directory}/labelled_minimap_alignment.paf"
labelled_minimap_alignment_new_paf="${output_directory}/labelled_minimap_alignment_new.paf"
flanks_bed="${output_directory}/flanks.bed"
flanks_fasta="${output_directory}/flanks.fasta"
flanks_fai="${output_directory}flanks.fasta.fai"
blast_db="${output_directory}/blast_sequences_db"
blast_output_dir="${output_directory}/blast_tables"



#Step 1: Run R script to go from input coordinates of genes to coordinates of 200bp anchors for each region

#echo "Running R to identify anchors"

#Rscript /g/data/te53/sj2852/blgrp/scripts/anchor_function.R -f "$input_bed_file" -o "$coordinate_anchors_bed"

if [[ ! -s "$coordinate_anchors_bed" ]]; then
    echo "Error: Anchors BED file not created: $coordinate_anchors_bed"
    exit 1
fi
echo "Anchors BED file generated successfully."

awk '
{
    col4 = $4
    for (i=5; i<=NF; i++) {
        col4 = col4 OFS $i
    }

    gsub(/[, ]/, "_", col4)

    print $1 "\t" $2 "\t" $3 "\t" col4
}' OFS=" " "$coordinate_anchors_bed" > "${coordinate_anchors_bed}.tmp"

mv "${coordinate_anchors_bed}.tmp" "$coordinate_anchors_bed"

sort -k1,1 -k2,2n "$coordinate_anchors_bed" -o "$coordinate_anchors_cleaned_bed"

if [[ ! -s "$coordinate_anchors_cleaned_bed" ]]; then
    echo "Error: Anchors BED file not created: $coordinate_anchors_cleaned_bed"
    exit 1
fi
echo "cleaned BED file generated successfully."

#Create a second BAM file that has the full sequences (anchor to anchor) for each region

awk '
{
    if ($4 != prev_gene) {
        if (NR > 1) {
            print prev_chr "\t" prev_start "\t" prev_end "\t" prev_gene
        }
        prev_chr = $1
        prev_start = $2
        prev_gene = $4
    }
    prev_end = $3
}
END {
    print prev_chr "\t" prev_start "\t" prev_end "\t" prev_gene
}' "$coordinate_anchors_cleaned_bed" > "$full_region_bed"

if [[ ! -s "$full_region_bed" ]]; then
    echo "Error: Full Region BED file not created: $full_region_bed"
    exit 1
fi
echo "Full Region BED file generated successfully."

#Step 2: Use Samtools to extract sequences of anchors

awk '{print $1":"$2"-"$3}' "$coordinate_anchors_cleaned_bed" | while read region; do samtools faidx "$reference_fasta" "$region"; done > "$blgrp_anchors_fasta"

if [[ ! -s "$blgrp_anchors_fasta" ]]; then
    echo "Error: Anchor FASTA file not created: $blgrp_anchors_fasta"
    exit 1
fi
echo "Anchor FASTA file generated successfully."

#Step 3.1: Use Minimap2 to align anchor sequences to chm13

minimap2 --cs "$reference_fasta" "$blgrp_anchors_fasta" > "$minimap_anchor_alignment_initial_paf"


if [[ ! -s "$minimap_anchor_alignment_initial_paf" ]]; then
    echo "Error: Initial Minimap file not created: $minimap_anchor_alignment_initial_paf"
    exit 1
fi
echo "Initial Minimap file generated successfully."

#Step 3.2: Use BLAST to align anchor sequences to chm13

makeblastdb -in "$reference_fasta" -dbtype nucl -out $blast_reference_db

blastn -query "$blgrp_anchors_fasta" -db "$blast_reference_db" -outfmt 6 -out "$blast_anchor_alignment_initial_tbl" 

awk '
{
    count[$1]++
    if (count[$1] <= 3) {
        print $0
    }
}' "$blast_anchor_alignment_initial_tbl" > "$anchor_top3_tbl"

awk '($3 ==100) || !($3 < 90 || $4 < 199 || $5 > 10)' "$anchor_top3_tbl" > "$anchor_filtered_tbl"

if [[ ! -s "$anchor_filtered_tbl" ]]; then
    echo "Error: Filtered Anchor file not created: $anchor_filtered_tbl"
    exit 1
fi
echo "Filtered anchor file generated successfully."


#Step 4: Label Genes 

awk '
BEGIN { OFS="\t" } 
{
    if (!seen[$1]++) {
        # Split the first column into chrom, gene_start, gene_end
        split($1, coord, "[:-]")
        chrom = coord[1]
        gene_start = coord[2]
        gene_end = coord[3]

        # Extract columns 9 and 10 (assumed to be $9 and $10)
        aln_start = $9
        aln_end = $10

        # Print chrom, gene_start, gene_end, aln_start, aln_end
        print chrom, gene_start, gene_end, aln_start, aln_end
    }
}' "$minimap_anchor_alignment_initial_paf" > "$labelled_minimap_alignment_paf"

awk '
BEGIN {
    while (getline < "'$coordinate_anchors_bed'") {
        bed_key = $1 ":" $2 "-" $3
        bed_genes[bed_key] = $4
    }
    close("'$coordinate_anchors_bed'")
}
{
    paf_key = $1 ":" $2 "-" $3
    if (paf_key in bed_genes) {
        print $0 "\t" bed_genes[paf_key]
    } else {
        print $0 "\t" "NA"
    }
}' "$labelled_minimap_alignment_paf" > "$labelled_minimap_alignment_new_paf"

if [[ ! -s "$labelled_minimap_alignment_new_paf" ]]; then
    echo "Error: labelled minimap file not created: $labelled_minimap_alignment_new_paf"
    exit 1
fi
echo "labelled minimap file generated successfully."

#Step 5: Add 300kb to either side of all genes

temp_file="/g/data/te53/sj2852/blgrp/tmp/full_coordinate_test/take_2/full_region.bed"

awk -v OFS="\t" -v prefix=300000 '
    {
        chrom = $1
        start = $2 - prefix
        end = $3 + prefix
        gene = $4

        if (start < 0) start = 0

        print chrom, start, end, gene
    }
' "$temp_file" > "$flanks_bed" 

if [[ ! -s "$flanks_bed" ]]; then
    echo "Error: Flanks BED file not created: $flanks_bed"
    exit 1
fi
echo "flanks BED file generated successfully."

rm "$temp_file"

awk '{print $1":"$2"-"$3}' "$flanks_bed" | while read region; do samtools faidx "$reference_fasta" "$region"; done > "$flanks_fasta"

if [[ ! -s "$flanks_fasta" ]]; then
    echo "Error: Flanks FASTA file not created: $flanks_fasta"
    exit 1
fi
echo "flanks FASTA file generated successfully."

#Step 6: Run Nucleotide BLAST on new coordinates against themselves

samtools faidx $flanks_fasta

for i in $(cut -f1 "$flanks_fai"); do
    
    samtools faidx "$flanks_fasta" "$i" > "$blast_output_dir/$i.fasta"

    makeblastdb -in "$blast_output_dir/$i.fasta" -dbtype nucl
    
    # Run BLAST and save the output to a .tbl file
    blastn -query "$blast_output_dir/$i.fasta" -db "$blast_output_dir/$i.fasta"  -outfmt 6 -out "$blast_output_dir/$i.self.alignment.tbl"

    rm "$blast_output_dir/$i.fasta.n*"
done

find "$blast_output_dir" -type f ! -name "*.fasta" ! -name "*.self.alignment.tbl" -exec rm {} +

#Step 7: Identify Unique Anchor Coordinates (update so it can be run)

Rscript /g/data/te53/sj2852/blgrp/scripts/anchor_scripts/anchor_idneitication_HP.R

#Step 8: Organise Output File

(head -n 1 /g/data/te53/sj2852/blgrp/tmp/rerunning_new_coords/NEW_COORDS/new_blgrp_coords.txt && tail -n +2 /g/data/te53/sj2852/blgrp/tmp/rerunning_new_coords/NEW_COORDS/new_blgrp_coords.txt | sort -k1,1V) > "/g/data/te53/sj2852/blgrp/tmp/rerunning_new_coords/NEW_COORDS/new_blgrp_coords_sorted.txt"

awk 'NR>1 {print $2, $5, $6; print $2, $7, $8}' /g/data/te53/sj2852/blgrp/tmp/rerunning_new_coords/NEW_COORDS/new_blgrp_coords_sorted.txt > "/g/data/te53/sj2852/blgrp/tmp/rerunning_new_coords/NEW_COORDS/new_blgrp_coords_sorted_transformed.txt"


#Get FASTA sequence for new bed reads

awk '{print $1":"$2"-"$3}' "$coordinate_anchors_cleaned_bed" | while read region; do samtools faidx "$reference_fasta" "$region"; done > "$blgrp_anchors_fasta"

#Run BLAST to confirm whole region is unique (confirm via top 3 alignments)

blastn -query /g/data/te53/sj2852/blgrp/tmp/rerunning_new_coords/NEW_COORDS/anchor_coords.fasta -db /g/data/te53/sj2852/blgrp/tmp/rerunning_new_coords/NEW_COORDS/blast/blast_reference_db -outfmt 6 -out /g/data/te53/sj2852/blgrp/tmp/rerunning_new_coords/NEW_COORDS/anchor_blast.tbl

awk '
> {
>     count[$1]++
>     if (count[$1] <= 3) {
>         print $0
>     }
> }' /g/data/te53/sj2852/blgrp/tmp/rerunning_new_coords/NEW_COORDS/anchor_blast.tbl > /g/data/te53/sj2852/blgrp/tmp/rerunning_new_coords/NEW_COORDS/anchor_blast_top3.tbl

#IF UNIQUE: Extract 200bp closest to gene body for each region

awk -v OFS="\t" '
    {
        # Output first row
        print $2, $6 - 200, $6
        
        # Output second row
        print $2, $7, $7 + 200
    }
' "/g/data/te53/sj2852/blgrp/tmp/rerunning_new_coords/NEW_COORDS/new_blgrp_coords_sorted.txt" > "/g/data/te53/sj2852/blgrp/tmp/rerunning_new_coords/NEW_COORDS/200bp_anchors/200bp_anchors.txt"

#Get FASTA of 200bp regions and BLAST again

awk '{print $1":"$2"-"$3}' "/g/data/te53/sj2852/blgrp/tmp/rerunning_new_coords/NEW_COORDS/200bp_anchors/200bp_anchors.txt" | while read region; do samtools faidx "/g/data/te53/sj2852/blgrp/referenceresource/chm13-t2t.ebv.phix.chrQ.xy.fa" "$region"; done > "/g/data/te53/sj2852/blgrp/tmp/rerunning_new_coords/NEW_COORDS/200bp_anchors/200bp_anchors.fasta"

blastn -query /g/data/te53/sj2852/blgrp/tmp/rerunning_new_coords/NEW_COORDS/200bp_anchors/200bp_anchors.fasta \
       -db /g/data/te53/sj2852/blgrp/tmp/rerunning_new_coords/NEW_COORDS/blast/blast_reference_db \
       -outfmt 6 \
       -out /g/data/te53/sj2852/blgrp/tmp/rerunning_new_coords/NEW_COORDS/200bp_anchors/200bp_anchor_blast.tbl

awk '
> {
>     count[$1]++
>     if (count[$1] <= 3) {
>         print $0
>     }
> }' /g/data/te53/sj2852/blgrp/tmp/rerunning_new_coords/NEW_COORDS/200bp_anchors/200bp_anchor_blast.tbl > /g/data/te53/sj2852/blgrp/tmp/rerunning_new_coords/NEW_COORDS/200bp_anchors/200bp_anchor_top3.tbl

#Generate new blgrp coordinates
