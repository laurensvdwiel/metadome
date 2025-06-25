import argparse

import os
import glob

import pandas as pd
from tqdm import tqdm
import gzip
from multiprocessing import Pool
from functools import partial


def parse_arguments():
    parser = argparse.ArgumentParser(description="Merge IDMapper, Metaposition and Genomic Coordinates")
    parser.add_argument("--input_csv", required=True)
    return parser.parse_args()

def main():
    
    args = parse_arguments()
    
    final = pd.read_csv(args.input_csv)

    
    #  Collect stats about this release
    # Calculate metrics
    print("Calculating stats of this release...")
    
    # Keep only needed columns
    final = final[[
    'chr', 'hg38', 'REF', 'ALT',
    'GeneSymbol', 'ENSEMBL_TR', 'uniprot_ac',
    'PFAM_ID', 'range_id', 'PFAM_consensus_pos', 'GencodeBasic' ]]

    metrics_current_release = {
        "Protein-coding genes": final['GeneSymbol'].nunique(),
        "Unique ENSEMBL protein coding transcripts": final['ENSEMBL_TR'].nunique(),
        "Unique SwissProt": final['uniprot_ac'].nunique(),
        "Unique PFAM": final['PFAM_ID'].nunique(),
        "Chromosome to protein position mappings": len(final),
        "Unique chromosome positions": final[['chr', 'hg38', 'REF', 'ALT']].drop_duplicates().shape[0],
        "Chromosome to protein position mappings with PFAM": final[final['PFAM_ID'].notna()].shape[0],
        "Pfam protein domain regions": final[['ENSEMBL_TR', 'range_id']].drop_duplicates().shape[0],
        "Unique ENSEMBL sequences with at least one Pfam domain annotated":
            final.loc[final['PFAM_ID'].notna(), 'ENSEMBL_TR'].nunique(),
        "Unique ENSEMBL sequences with at least two Pfam domain annotated":
            final[final['PFAM_ID'].notna()]
                .drop_duplicates(subset=['ENSEMBL_TR', 'range_id'])
                .groupby('ENSEMBL_TR')
                .filter(lambda x: len(x) >= 2)['ENSEMBL_TR']
                .nunique(),
        "Unique SwissProt sequences with at least one Pfam domain annotated":
            final.loc[final['PFAM_ID'].notna(), 'uniprot_ac'].nunique(),
        "Unique SwissProt sequences with at least two Pfam domain annotated":
            final[final['PFAM_ID'].notna()]
                .drop_duplicates(subset=['uniprot_ac', 'range_id'])
                .groupby('uniprot_ac')
                .filter(lambda x: len(x) >= 2)['uniprot_ac']
                .nunique(),
        "Average number of homologs per Pfam": final[final['PFAM_ID'].notna()]
            .drop_duplicates(subset=['PFAM_ID', 'ENSEMBL_TR', 'range_id'])
            .groupby('PFAM_ID').size().mean(),
        "Average length of Pfam domains": final.dropna(subset=['PFAM_consensus_pos'])
            .groupby(['ENSEMBL_TR', 'range_id'])['PFAM_consensus_pos'].max().mean()
    }

    # Subset for basic gencode to compare with initial metadome
    final_basic = final[final['GencodeBasic']==True]

    metrics_current_release.update({
        "Protein-coding genes (GencodeBasic)": final_basic['GeneSymbol'].nunique(),
        "Unique ENSEMBL protein coding transcripts (GencodeBasic)": final_basic['ENSEMBL_TR'].nunique(),
        "Unique SwissProt (GencodeBasic)": final_basic['uniprot_ac'].nunique(),
        "Unique PFAM (GencodeBasic)": final_basic['PFAM_ID'].nunique(),
        "Chromosome to protein position mappings (GencodeBasic)": len(final_basic),
        "Unique chromosome positions (GencodeBasic)":
            final_basic[['chr', 'hg38', 'REF', 'ALT']].drop_duplicates().shape[0],
        "Chromosome to protein position mappings with PFAM (GencodeBasic)":
            final_basic[final_basic['PFAM_ID'].notna()].shape[0],
        "Pfam protein domain regions (GencodeBasic)":
            final_basic[['ENSEMBL_TR', 'range_id']].drop_duplicates().shape[0],
        "Unique ENSEMBL sequences with at least one Pfam domain annotated (GencodeBasic)":
            final_basic.loc[final_basic['PFAM_ID'].notna(), 'ENSEMBL_TR'].nunique(),
        "Unique ENSEMBL sequences with at least two Pfam domain annotated (GencodeBasic)":
            final_basic[final_basic['PFAM_ID'].notna()]
                .drop_duplicates(subset=['ENSEMBL_TR', 'range_id'])
                .groupby('ENSEMBL_TR')
                .filter(lambda x: len(x) >= 2)['ENSEMBL_TR']
                .nunique(),
        "Unique SwissProt sequences with at least one Pfam domain annotated (GencodeBasic)":
            final_basic.loc[final_basic['PFAM_ID'].notna(), 'uniprot_ac'].nunique(),
        "Unique SwissProt sequences with at least two Pfam domain annotated (GencodeBasic)":
            final_basic[final_basic['PFAM_ID'].notna()]
                .drop_duplicates(subset=['uniprot_ac', 'range_id'])
                .groupby('uniprot_ac')
                .filter(lambda x: len(x) >= 2)['uniprot_ac']
                .nunique(),
        "Average number of homologs per Pfam (GencodeBasic)":
            final_basic[final_basic['PFAM_ID'].notna()]
                .drop_duplicates(subset=['PFAM_ID', 'ENSEMBL_TR', 'range_id'])
                .groupby('PFAM_ID').size().mean(),
        "Average length of Pfam domains (GencodeBasic)":
            final_basic.dropna(subset=['PFAM_consensus_pos'])
                .groupby(['ENSEMBL_TR', 'range_id'])['PFAM_consensus_pos'].max().mean()
    })


    pd.DataFrame(list(metrics_current_release.items()), columns=['Stat', 'Counts']).to_csv(args.input_csv + "stats", index=False)

if __name__ == "__main__":

    main()

