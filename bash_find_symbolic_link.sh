# This script to find and create symbolic links for certain type of files.
# Created by Hyungtaek Jung

#!/bin/bash

# Usage: bash bash_find_symbolic_link.sh

# Step 2: Set the present work directory
CURRENT_FOLDER="/g/data/te53/sj2852/blgrp/data/KGPdata/kgp_all_fastq"

# Step 3: Set the file patterns to find
#FIND_PATTERN=(".flyasm.fasta" ".flyasm.fasta.gz" ".flyasm.fa" ".flyasm.fa.gz")
FIND_PATTERN=(".blgrp.fastq")

# Step 4: Set the output directory for symbolic links
OUTDIR=/g/data/te53/sj2852/blgrp/analyses/flye_assemblies/hapdup/fastq_files

# Ensure the output directory exists
mkdir -p "$OUTDIR"

# Function to create symbolic links
create_symlinks() {
    local pattern="$1"
    local found_files

    # Find files matching the pattern
    found_files=$(find "$CURRENT_FOLDER" -type f -name "*$pattern")

    for file in $found_files; do
        # Get the base name of the file
        local basename=$(basename "$file")
        # Create the symbolic link in the output directory
        ln -s "$file" "$OUTDIR/$basename"
        echo "Created symbolic link for $file in $OUTDIR/$basename"
    done
}

# Iterate over the find patterns and create symbolic links
for pattern in "${FIND_PATTERN[@]}"; do
    create_symlinks "$pattern"
done

echo "Symbolic link creation completed."
