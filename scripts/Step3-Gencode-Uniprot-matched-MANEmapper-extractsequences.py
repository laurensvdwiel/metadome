import pandas as pd
import numpy as np
import argparse
import gzip
from Bio import SeqIO

def parse_arguments():
    parser = argparse.ArgumentParser(description="Annotate and filter FASTA file using a file containing Gencode-Uniprot matches and add MANE annotations.")
    parser.add_argument("--mapper_csv", required=True, help="Path to the mapping CSV file.")
    parser.add_argument("--protein_fasta", required=True, help="Path to the protein FASTA file.")
    parser.add_argument("--transcript_fasta", required=True, help="Path to the transcript FASTA file.")
    parser.add_argument("--output_csv", required=True, help="Output path for the updated CSV.")
    parser.add_argument("--output_protein_fa", required=True, help="Output path for the filtered protein FASTA.")
    
    return parser.parse_args()


def parse_fasta(fasta_file, protein_ids):
    if fasta_file.endswith('.gz'):
        with gzip.open(fasta_file, 'rt') as handle:  # 'rt' mode for text reading from gzip
            filtered_protein_seqs = (record for record in SeqIO.parse(handle, "fasta") 
                                     if record.id in protein_ids)
            return list(filtered_protein_seqs)
    else:
        filtered_protein_seqs = (record for record in SeqIO.parse(fasta_file, "fasta") 
                                if record.id in protein_ids)
        return filtered_protein_seqs


def main():
    
    args = parse_arguments()
    
    # Load and process the mapper file
    mappers = pd.read_csv(args.mapper_csv)
    mappers['Ensembl_PRT'] = mappers['qseqid'].str.split("|", expand=True)[0]
    mappers['Gene'] = mappers['qseqid'].str.split("|", expand=True)[7]
    
    # Save the updated mapper file
    mappers.to_csv(args.output_csv, index=False)

    # Filter and write protein FASTA
    protein_ids = set(mappers.qseqid)
    filtered_protein_seqs = parse_fasta(args.protein_fasta, protein_ids)
    SeqIO.write(filtered_protein_seqs, args.output_protein_fa, "fasta")

if __name__ == "__main__":
    main()
