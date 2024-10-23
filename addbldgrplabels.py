import pandas as pd
import argparse
import gzip
from Bio import SeqIO
from Bio.Seq import Seq


def load_paf(file_path):
    """
    Load PAF file data into a pandas DataFrame.
    """
    columns = [
        'query_name', 'query_length', 'query_start', 'query_end', 'relative_strand',
        'target_name', 'target_length', 'target_start_on_strand', 'target_end_on_strand',
        'residue_matches'
    ]
    return pd.read_csv(file_path, sep='\t', header=None, names=columns, usecols=range(len(columns)))


def read_fasta_headers(fasta_file):
    """
    Read headers from a FASTA file. Supports both gzipped and plain text files.
    """
    headers = []
    open_func = gzip.open if fasta_file.endswith('.gz') else open

    with open_func(fasta_file, 'rt') as file:
        for line in file:
            if line.startswith('>'):
                header = line[1:].strip().split()[0]
                headers.append(header)
    return headers


def write_modified_fasta(fasta_file, merged_data, output_file):
    """
    Write a modified FASTA file with the label appended as #{label} in the header.
    Sequences are reversed and complemented if the strand is negative, and the 'N:flipped' tag is added.
    """
    open_func = gzip.open if fasta_file.endswith('.gz') else open

    with open_func(fasta_file, 'rt') as file_in, open(output_file, 'w') as file_out:
        fasta_sequences = SeqIO.parse(file_in, 'fasta')

        for fasta in fasta_sequences:
            # Get sequence id and corresponding metadata from merged_data dataframe
            query_name = fasta.id
            sequence = str(fasta.seq)

            # Fetch the corresponding row from the dataframe
            row = merged_data[merged_data['query_name'] == query_name].iloc[0]
            label = row['label']
            strand = row['relative_strand']

            # Modify the header with #{label}
            new_header = f">{fasta.id}#{label}"

            # Check if the sequence needs to be flipped
            if strand == '-':
                sequence = str(Seq(sequence).reverse_complement())  # Flip the sequence
                new_header += " N:flipped"  # Add the N: tag

            # Write the modified header and sequence to the new file
            file_out.write(new_header + "\n")
            file_out.write(sequence + "\n")


def main():
    # Set up argparse for command-line arguments
    parser = argparse.ArgumentParser(description="Process PAF and FASTA files to modify headers and sequences.")
    parser.add_argument('--paf', required=True, help="Input PAF file")
    parser.add_argument('--genesbed', required=True, help="Input BED file with gene information")
    parser.add_argument('--inputfasta', required=True, help="Input FASTA file")
    parser.add_argument('--outputfasta', required=True, help="Output FASTA file")
    args = parser.parse_args()

    # Load PAF data
    paf_data = load_paf(args.paf)

    # First DataFrame: Calculate summed residue matches percentage per query_name and target_name
    grouped_paf_data = paf_data.groupby(['query_name', 'target_name']).agg({
        'residue_matches': 'sum',
        'query_length': 'first',
        'target_length': 'first'
    }).reset_index()

    grouped_paf_data['residue_match_percentage'] = (grouped_paf_data['residue_matches'] / grouped_paf_data['query_length']) * 100

    # Keeping the row with the highest residue_match_percentage for each query_name
    best_matches = grouped_paf_data.loc[grouped_paf_data.groupby('query_name')['residue_match_percentage'].idxmax()]

    # Second DataFrame: Get max residue matches for each target_name and query_name pair
    max_residue_matches = paf_data.loc[paf_data.groupby(['target_name', 'query_name'])['residue_matches'].idxmax()]
    max_residue_matches = max_residue_matches[['target_name', 'query_name', 'query_start', 'query_end', 'query_length', 'target_length', 'relative_strand']]

    # Merge the two DataFrames
    merged_data = pd.merge(
        best_matches[['target_name', 'query_name', 'residue_match_percentage']],
        max_residue_matches,
        on=['target_name', 'query_name']
    )

    # Round percentages and target_start_on_strand for cleaner output
    merged_data['residue_match_percentage'] = round(merged_data['residue_match_percentage'], 2)

    # Load gene information from BED file and merge with the PAF data
    genes_bed = pd.read_csv(args.genesbed, sep='\t')
    genes_bed['target_name'] = genes_bed['region']
    genes_bed['gene'] = genes_bed['gene'].str.replace(r'-[53]$', '', regex=True)
    genes_bed = genes_bed[['target_name', 'gene']]
    genes_bed = genes_bed.drop_duplicates()

    merged_data = pd.merge(merged_data, genes_bed, on=['target_name'])
    merged_data['label'] = merged_data['gene']

    # Read the FASTA headers and prepare final merged DataFrame
    fasta_headers = read_fasta_headers(args.inputfasta)
    fasta_df = pd.DataFrame(fasta_headers, columns=['query_name'])

    final_merged_data = pd.merge(fasta_df, merged_data, on='query_name', how='left')

    # Assign geneUn label to contigs present in FASTA but not in the PAF results
    final_merged_data['label'] = final_merged_data['label'].fillna(final_merged_data['target_name'].apply(lambda x: 'geneUn'))

    # Write the modified FASTA file
    write_modified_fasta(args.inputfasta, final_merged_data, args.outputfasta)


if __name__ == "__main__":
    main()
