import gffutils
import pandas as pd
import pyfaidx
from Bio.Data import CodonTable
from Bio.Seq import Seq
from scipy.special import softmax
import argparse
import os
import pickle
import multiprocessing
from functools import partial

def parse_arguments():
    parser = argparse.ArgumentParser(description="Computes all possible single nucleotide variants for a given set of transcripts.")
    parser.add_argument('--gtf', required=True, help='Path to GTF file (.gz)')
    parser.add_argument('--genome_fa', required=True, help='Path to genome fasta file')
    parser.add_argument('--output_dir', required=True, help='Directory to save annotated variants')
    parser.add_argument('--db_path', default='data/intermediate/gencode.annotation.gtf.db', help='Path to GTF database file')
    parser.add_argument('--n_cores', type=int, default=None, help='Number of CPU cores to use')
    parser.add_argument('--transcript_csv', default=None, help='Optional CSV file with column transcript_id to filter transcripts')
    return parser.parse_args()



def extract_transcripts_by_gene(db, fasta_file):
    """
    Extracts transcript information for protein-coding genes from a genomic database (gffutils-parsed GTF) and a FASTA file.
    """
    # Load the FASTA file
    genome = pyfaidx.Fasta(fasta_file)

    # Dictionary to store transcript information
    transcripts_dict = {}

    # Get all genes from the database
    for gene in db.features_of_type("gene"):

        if gene.attributes.get('gene_type')[0] == "protein_coding":

            # Get all transcripts for this gene
            transcripts = db.children(gene, featuretype='transcript')

            for transcript in transcripts:
                transcript_id = transcript.id
                chromosome = transcript.seqid
                strand = transcript.strand
                
                # Get all tags 
                tags = transcript.attributes.get("tag", [])
                # Gencode basic
                is_basic = "basic" in tags
                # Handling MANE, expecting a tag with MANE_Select MANE_Plus_Clinical or none
                mane_tag = next((tag for tag in tags if tag.startswith("MANE")), "-")
                
                # Get all exons for this transcript and sort them
                exons = list(db.children(transcript, featuretype='exon'))
                exons.sort(key=lambda x: x.start)

                # Holder for the exons 
                exon_regions = []

                # Extract and concatenate exon sequences to get the full transcript
                transcript_sequence = ""
                for exon in exons:
                    seq = genome[chromosome][exon.start - 1:exon.end].seq
                    transcript_sequence += seq
                    exon_regions.append({
                        'start': exon.start,
                        'end': exon.end,
                        'exon_number': exon.attributes.get('exon_number', ['-'])[0],
                        'tags': exon.attributes.get('tag', []),  # can be empty
                    })

                # Reverse complement if on negative strand
                if strand == '-':
                    transcript_sequence = str(Seq(transcript_sequence).reverse_complement())

                # Also get CDS features to determine coding sequences
                cds_features = list(db.children(transcript, featuretype='CDS'))
                cds_features.sort(key=lambda x: x.start)
                cds_regions = [(cds.start, cds.end) for cds in cds_features]

                # Extract coding sequence if present
                coding_sequence = ""
                if cds_features:
                    for cds in cds_features:
                        seq = genome[chromosome][cds.start - 1:cds.end].seq
                        coding_sequence += seq

                    if strand == '-':
                        coding_sequence = str(Seq(coding_sequence).reverse_complement())

                # Create dictionary entry
                transcripts_dict[transcript_id] = {
                    'transcriptID': transcript_id,
                    'geneID': gene.id,
                    'chromosome': chromosome,
                    'transcript_length': len(transcript_sequence),
                    'coding_length': len(coding_sequence) if coding_sequence else 0,
                    'strand': strand,
                    'exon_regions': exon_regions,
                    'cds_regions': cds_regions if cds_features else [],
                    'transcript_sequence': transcript_sequence,
                    'coding_sequence': coding_sequence if coding_sequence else "",
                    'is_coding': bool(cds_features),
                    'is_basic': is_basic,
                    'mane_tag': mane_tag
                }

    return transcripts_dict


def calculate_genomic_position(transcript_pos, strand, cds_regions):
    """
    Calculate the genomic position corresponding to a given transcript position.
    This function maps a position within a transcript to its corresponding genomic
    position based on the coding sequence (CDS) regions and strand orientation.
    """
    
    # Track cumulative length and find the positions
    cumulative_length = 0

    cds_regions_oredered = cds_regions if strand == '+' else cds_regions[::-1]
    
    for region in cds_regions_oredered:

        region_length = region[1] - region[0] + 1

        if cumulative_length + region_length > transcript_pos:
            # Position is in this region
            if strand == '+':
                # For positive strand, add transcript position to region start
                return region[0] + (transcript_pos - cumulative_length)
            else:
                # For negative strand, subtract from region end
                return region[1] - (transcript_pos - cumulative_length)

        cumulative_length += region_length

    # If position not found (should not happen)
    raise ValueError("Transcript position out of bounds")


def generate_codon_table(transcript_info):
    """
   This function takes a transcript's coding sequence and generates a table with all possible codons 
    """
    transcript_seq = transcript_info['coding_sequence']

    variants = []

    # Use standard genetic code from Biopython
    genetic_code = CodonTable.standard_dna_table.forward_table
    stop_codons = CodonTable.standard_dna_table.stop_codons

    # Strand information
    strand = transcript_info['strand']
    cds_regions = transcript_info['cds_regions']

    # Iterate through each position
    for pos in range(len(transcript_seq)):
        
        ref_base = transcript_seq[pos]

        # Codon and amino acid calculations
        codon_index = pos // 3
        codon_start = codon_index * 3
        codon_end = codon_start + 3

        ref_codon = transcript_seq[codon_start:codon_end]
        
        # Handle amino acid translation
        # Explicitly check for stop codons
        if ref_codon in stop_codons:
            ref_aa = 'X'
        else:
            ref_aa = genetic_code.get(ref_codon, 'X') if len(ref_codon) == 3 else 'X'

        # Calculate genomic position
        genomic_pos = calculate_genomic_position(pos, strand, cds_regions)

        # Determine which exon the genomic position belongs to
        exon_number = "-"
        for exon in transcript_info['exon_regions']:
            if exon['start'] <= genomic_pos <= exon['end']:
                exon_number = exon['exon_number']
                break

        variants.append({
            'transcriptID': transcript_info.get('transcriptID', 'None'),
            'chr': transcript_info.get('chromosome', 'None'),
            'hg38': genomic_pos,
            'REF': ref_base, #if strand == '+' else str(Seq(ref_base).reverse_complement()),
            'codon' : ref_codon,
            'cDNA_position': pos + 1,
            'codon_base_pair_position': pos - codon_start,
            'strand': strand,
            'RefAA': ref_aa,
            'AAindex': codon_index + 1,
            'exon_number': exon_number,
            'GencodeBasic': transcript_info.get('is_basic', False),
            'MANE': transcript_info.get('mane_tag', '-')
        })
    return pd.DataFrame(variants)


def process_transcript(trn, transcripts_dict, output_dir):
    output_path = os.path.join(output_dir, f"{trn}.tsv")

    if os.path.exists(output_path):
        print(f"Skipping {trn}, as already annotated ")
        return

    print(f"Computing possible variants for {trn}")
    transcript_info = transcripts_dict.get(trn)

    if transcript_info and transcript_info.get('coding_sequence'):
        variants = generate_codon_table(transcript_info)
        variants.to_csv(output_path, sep='\t', index=False)


def main():
    args = parse_arguments()

    os.makedirs(args.output_dir, exist_ok=True)

    if not os.path.exists(args.db_path):
        print("Creating transcript database")
        db = gffutils.create_db(args.gtf, dbfn=args.db_path, force=True,
                                merge_strategy='merge',
                                sort_attribute_values=False,
                                keep_order=False,
                                verbose=True,
                                disable_infer_transcripts=True,
                                disable_infer_genes=True)
    else:
        print("Loading transcript database")
        db = gffutils.FeatureDB(args.db_path, keep_order=True)

    # Load or create transcript dictionary
    transcript_cache = "data/intermediate/gencode.transcripts.cache.pkl"
    if os.path.exists(transcript_cache):
        print("Loading cached transcript dictionary...")
        with open(transcript_cache, "rb") as f:
            transcripts_dict = pickle.load(f)
    else:
        print("Extracting transcripts...")
        transcripts_dict = extract_transcripts_by_gene(db, args.genome_fa)
        with open(transcript_cache, "wb") as f:
            pickle.dump(transcripts_dict, f)

    # Filter transcripts if CSV provided
    if args.transcript_csv:
        trn_df = pd.read_csv(args.transcript_csv)        
        trn_ids = set(trn_df['transcript_id'].str.split('|', expand=True)[0])
        transcripts_dict = {k: v for k, v in transcripts_dict.items() if k in trn_ids}
        print(f"Filtered to {len(transcripts_dict)} transcripts based on input list")

    # Determine cores
    with multiprocessing.Pool(args.n_cores) as pool:
        pool.map(partial(process_transcript, transcripts_dict=transcripts_dict, output_dir=args.output_dir),
                 transcripts_dict.keys())


if __name__ == "__main__":
    main()