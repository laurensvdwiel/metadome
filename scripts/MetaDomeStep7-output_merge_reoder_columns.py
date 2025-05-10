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
    parser.add_argument("--uniprot_name", required=True)
    parser.add_argument("--pfamscan_output", required=True)
    parser.add_argument("--genome_build", required=True)
    parser.add_argument("--PFAM_version", required=True)
    parser.add_argument("--PFAM_interpro", required=False, default=None)
    parser.add_argument("--source", required=True, default=None)
    parser.add_argument("--GENCODE_version", required=True, default=None)
    parser.add_argument("--output", required=True)
    parser.add_argument('--n_cores', type=int, default=None, help='Number of CPU cores to use')
    return parser.parse_args()

def load_pfam_to_interpro(path):
    return pd.read_csv(path, 
                       usecols=['Accession', 'Integrated Into'], 
                       dtype=str).rename(columns={'Accession': 'PFAM_ID', 'Integrated Into': 'interpro_id'})

def load_idmapper(path):
    return pd.read_csv(path, 
                       usecols=['qseqid', 'sseqid'], 
                       dtype=str).rename(columns={'qseqid': 'ENSEMBL', 'sseqid': 'uniprot_ac'})

def load_metaposition(path):
    return pd.read_csv(path, sep='\t').rename(columns={
        'transcript_id': 'ENSEMBL',
        'domain_id': 'PFAM_ID',
        'range_id': 'range_id',
        'seq_position': 'uniprot_pos',
        'sequence_aa': 'uniprot_AA',
        'position_in_domain_consensus': 'PFAM_consensus_pos'
    })


def load_gencode_refseq(path):
    with gzip.open(path, 'rt') as f:
        return pd.read_csv(f, sep="\t", names=['ENSEMBL_TR', "RefSeq", "RefSeq_prot"],
                           usecols=[0, 1], dtype=str)


def parse_pfamscan(pfamscan_path):
    """
    Parse pfamscan output from disk into a DataFrame of unique PFAM_ID + name pairs.
    """
    if not os.path.isfile(pfamscan_path):
        raise FileNotFoundError(f"{pfamscan_path} does not exist")
    with open(pfamscan_path, 'r') as f:
        lines = f.read().splitlines()

    data_lines = [l for l in lines if l and not l.startswith('#')]

    cols = [
        'seq_id','align_start','align_end','env_start','env_end',
        'hmm_acc','hmm_name','type','hmm_start','hmm_end',
        'hmm_length','bit_score','e_value','significance','clan'
    ]

    parsed = []
    for line in data_lines:
        parts = line.split()
        if len(parts) != len(cols):
            continue
        parsed.append(dict(zip(cols, parts)))

    df = pd.DataFrame(parsed)
    df = (
        df.assign(
            PFAM_ID = df['hmm_acc'].str.split('.').str[0],
            name    = df['hmm_name']
        )
        [['PFAM_ID','name']]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    return df


def split_ensembl_fields(df):
    split_fields = df['ENSEMBL'].str.split('|', expand=True)
    df['ENSEMBL_TR'] = split_fields[1]  
    df['gencode_gene_id'] = split_fields[2]  
    df['gencode_translation_name'] = split_fields[5]  
    df['havana_gene_id'] = split_fields[3]  
    df['havana_translation_id'] = split_fields[4] 
    df['sequence_length'] = split_fields[7]  
    df['gene_name'] = split_fields[6]  
    
    return df
        

def _load_single_coordinate(folder, tr):
    file_path = os.path.join(folder, f"{tr}.tsv")
    if os.path.isfile(file_path):
        df = pd.read_csv(file_path, sep='\t')
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

    # Add uniprot names
    uniprot_name = pd.read_csv(args.uniprot_name, sep="|", names=['sp', "uniprot_ac", "uniprot_name"])
    idmapper = pd.merge(idmapper, uniprot_name, on='uniprot_ac', how='inner')
    del uniprot_name
    
    # Split ENSEMBL fields
    idmapper = split_ensembl_fields(idmapper)
    metaposition = load_metaposition(args.metaposition)
    # Add positions
    metaposition['uniprot_pos'] = metaposition['uniprot_pos'] - 1
    metaposition['uniprot_start'] = metaposition.groupby(['ENSEMBL', 'PFAM_ID', 'range_id'])['uniprot_pos'].transform('min')
    metaposition['uniprot_stop'] = metaposition.groupby(['ENSEMBL', 'PFAM_ID', 'range_id'])['uniprot_pos'].transform('max')
    metaposition['domain_length'] = metaposition['uniprot_stop'] - metaposition['uniprot_start']
    
    # Add pfam domain names
    pfam_names = parse_pfamscan(args.pfamscan_output)
    metaposition = pd.merge(metaposition, pfam_names, on='PFAM_ID', how='right')
    del pfam_names
    
    # Add pfam-interpro mappings
    interpro = load_pfam_to_interpro(args.PFAM_interpro)
    metaposition = pd.merge(metaposition, interpro, on="PFAM_ID", how='left')
    
    # Add RefSeq
    refseq = load_gencode_refseq(args.refseq)
    merged = pd.merge(idmapper, refseq, on='ENSEMBL_TR', how="left")
    del idmapper, refseq
    
    # Load relevant genomic coordinate files
    ensembl_tr_set = merged['ENSEMBL_TR'].dropna().unique()
    genomic_chunks = load_genomic_coordinates(args.genomic_folder, ensembl_tr_set, args.n_cores)
    genomic = pd.concat(genomic_chunks, ignore_index=True)
    del genomic_chunks
    
    final = pd.merge(
        merged,
        genomic,
        on='ENSEMBL_TR',
        how='right'
    )
    del merged, genomic
    
    final['uniprot_pos'] = final['uniprot_pos'] -1 
    
    final = final.merge(metaposition, on=['ENSEMBL', 'uniprot_AA', 'uniprot_pos'], how='outer')
    del metaposition
    
    # Add overall info if a domain exists at all
    final['evaluated_interpro_domains'] = final.groupby('ENSEMBL_TR')['PFAM_ID'] \
                                          .transform(lambda x: x.any())
    

    final['genome_build'] = args.genome_build
    final['PFAM_version'] = args.PFAM_version
    final['GENCODE_version'] = args.GENCODE_version
    final['source'] = args.source
        
    # # Reorder columns
    final = final[[
        'chr', 'hg38', 'REF', 'strand',
        'gene_name', 'gencode_gene_id', 'havana_gene_id', 'exon_number',
        'gencode_translation_name','ENSEMBL_TR','RefSeq', 'havana_translation_id',
        'cDNA_position', 'codon', 'codon_base_pair_position',
        'uniprot_ac', 'uniprot_pos', 'uniprot_AA', 'sequence_length',
        'evaluated_interpro_domains',
        'PFAM_ID','name', 'interpro_id', 'domain_length' ,'PFAM_consensus_pos', 'uniprot_start', 'uniprot_stop',
        'MANE','GencodeBasic',
        'genome_build', 'PFAM_version','GENCODE_version','source']].rename(
            columns={'ENSEMBL_TR': 'gencode_transcription_id',
                       'uniprot_ac': 'uniprot_name',
                       'uniprot_pos': 'amino_acid_position',
                       'chr': 'chromosome',
                       'hg38': 'chromosome_position',
                       'REF': 'base_pair',
                       'uniprot_AA': 'amino_acid_residue',
                       'PFAM_ID': 'ext_db_id'})
        
    int_cols = [
        'chromosome_position','exon_number','cDNA_position',
        'codon_base_pair_position','amino_acid_position',
        'domain_length','PFAM_consensus_pos','uniprot_start','uniprot_stop'
    ]
    final[int_cols] = final[int_cols].astype('Int32')
    
    final.sort_values(['chromosome','chromosome_position'], inplace=True)
    
    # Add duplicate columns
    final['uniprot_position'] = final['amino_acid_position']
    final['uniprot_residue'] = final['amino_acid_residue']
    
    final.to_csv(args.output, index=False)
   
if __name__ == "__main__":

    main()

