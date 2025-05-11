import os
import glob
import argparse
import pandas as pd
from tqdm import tqdm
import gzip
from multiprocessing import Pool
from functools import partial


def parse_arguments():
    parser = argparse.ArgumentParser(description="Merge IDMapper, Metaposition and Genomic Coordinates")
    parser.add_argument("--idmapper", required=True)
    parser.add_argument("--metaposition", required=True)
    parser.add_argument("--genomic_folder", required=True)
    parser.add_argument("--refseq", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument('--n_cores', type=int, default=None, help='Number of CPU cores to use')
    return parser.parse_args()

def load_idmapper(path):
    return pd.read_csv(path, usecols=['qseqid', 'sseqid'], dtype=str).rename(columns={
        'qseqid': 'ENSEMBL', 'sseqid': 'uniprot_ac'})

def load_metaposition(path):
    df = pd.read_csv(path, sep='\t', dtype={
        'transcript_id': 'category',
        'domain_id': 'category',
        'range_id': 'category',
        'seq_position': 'Int32',
        'sequence_aa': 'string',
        'position_in_domain_consensus': 'Int32'
    }, usecols=[
        'transcript_id', 'domain_id', 'range_id',
        'seq_position', 'sequence_aa', 'position_in_domain_consensus'
    ]).rename(columns={
        'transcript_id': 'ENSEMBL',
        'domain_id': 'PFAM_ID',
        'range_id': 'range_id',
        'seq_position': 'uniprot_pos',
        'sequence_aa': 'uniprot_AA',
        'position_in_domain_consensus': 'PFAM_consensus_pos'
    })
    df['ENSEMBL_TR']  = df['ENSEMBL'].str.split('|', expand=True)[1]
    return df


def load_gencode_refseq(path):
    with gzip.open(path, 'rt') as f:
        return pd.read_csv(f, sep="\t", names=['ENSEMBL_TR', "RefSeq", "RefSeq_prot"],
                           usecols=[0, 1], dtype=str)


def split_ensembl_fields(df):
    split_fields = df['ENSEMBL'].str.split('|', expand=True)
    df['ENSEMBL_TR'] = split_fields[1]  
    df['GeneSymbol'] = split_fields[6]  
    return df
        

def _load_single_coordinate(folder, tr):
    file_path = os.path.join(folder, f"{tr}.tsv")
    if os.path.isfile(file_path):
        df = pd.read_csv(file_path, sep='\t', dtype={
            'chr': 'category', 'hg38': 'Int32', 'REF': 'category', 'ALT': 'category',
            'transcriptID': 'string', 'RefAA': 'string', 'AAindex': 'Int32','strand': 'category',
            'MANE': 'category', 'GencodeBasic': 'category', 'exon_number': 'Int16'
        }, usecols=lambda c: c in {
            'chr', 'hg38', 'REF', 'ALT', 'transcriptID', 'RefAA', 'AAindex',
            'MANE', 'GencodeBasic', 'exon_number', 'strand'
        })
        return df.rename(columns={
            'transcriptID': 'ENSEMBL_TR',
            'RefAA': 'uniprot_AA',
            'AAindex': 'uniprot_pos'
        })
    return None


def load_genomic_coordinates(folder, ensembl_tr_set, n_cores=1):
    with Pool(n_cores) as pool:
        results = list(tqdm(pool.imap_unordered(partial(_load_single_coordinate, folder), ensembl_tr_set),
                            total=len(ensembl_tr_set), desc='Loading genomic coordinates'))
    return (df for df in results if df is not None)


def main():
    
    args = parse_arguments()

    # Load files
    idmapper = load_idmapper(args.idmapper)
    # Split ENSEMBL fields
    idmapper = split_ensembl_fields(idmapper)
    # Add RefSeq
    refseq = load_gencode_refseq(args.refseq)
    merged = pd.merge(idmapper, refseq, on='ENSEMBL_TR', how="left")
    del refseq
    
    print(merged.columns)
    ensembl_tr_set = merged['ENSEMBL_TR'].dropna().unique()
    genomic_chunks = load_genomic_coordinates(args.genomic_folder, ensembl_tr_set, args.n_cores)
    genomic = pd.concat(genomic_chunks, ignore_index=True)
    del genomic_chunks
    print(genomic.columns)

    final = pd.merge(
        genomic,
        merged,
        on=['ENSEMBL_TR'],
        how='left'
    )
    del merged, genomic
    print(final.columns)
    
    #  Add metaposition info
    metaposition = load_metaposition(args.metaposition)

    print(metaposition.columns)
    

    final = pd.merge(
        final,
        metaposition,
        on=['ENSEMBL_TR', 'uniprot_AA', 'uniprot_pos'],
        how='outer'
    )
 
    # Reorder columns
    final = final[[
        'chr', 'hg38', 'REF', 'ALT', 'strand', 
        'GeneSymbol', 'ENSEMBL_TR', 'RefSeq', 'exon_number', 'uniprot_ac', 
        'uniprot_pos', 'uniprot_AA', 'PFAM_ID', 'range_id','PFAM_consensus_pos','MANE','GencodeBasic'
        
    ]]

    # Cast to integer the numeric columns and save
    final[["hg38", "uniprot_pos", "PFAM_consensus_pos"]] = final[["hg38", "uniprot_pos", "PFAM_consensus_pos"]].astype("Int64")
    
    print("Saving merged file...")
    final.to_csv(args.output, index=False)
    
    #  Collect stats about this release
    # Calculate metrics
    print("Calculating stats of this release...")
    
    # Keep only needed columns
    final = final[[
    'chr', 'hg38', 'REF', 'ALT',
    'GeneSymbol', 'ENSEMBL_TR', 'uniprot_ac',
    'PFAM_ID', 'range_id', 'PFAM_consensus_pos', ]]

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
    final_basic = final[final['GencodeBasic']]

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


    pd.DataFrame(list(metrics_current_release.items()), columns=['Stat', 'Counts']).to_csv(args.output + "stats", index=False)

if __name__ == "__main__":

    main()

