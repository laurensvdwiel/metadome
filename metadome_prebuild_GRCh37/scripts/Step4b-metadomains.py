import argparse
import gzip
import tempfile
import pandas as pd
from Bio import SeqIO
import subprocess
import shutil
import os
from concurrent.futures import ProcessPoolExecutor
import tempfile

# Need to add this to prevent micromamba conflicts when running this script 
os.environ["MAMBA_ROOT_PREFIX"] = tempfile.mkdtemp()

def parse_arguments():
        parser = argparse.ArgumentParser(description="Extract metadomain information from PfamScan output and prepares files to compute metapositions.")
        parser.add_argument("--pfamscan", required=True)
        parser.add_argument("--fasta", required=True)
        parser.add_argument("--hmmdatabase", required=True)
        parser.add_argument("--output_dir", required=True)
        parser.add_argument("--padding", type=int, default=3, help="Padding for domain extraction (default=3)")
        parser.add_argument("--n_cores", type=int, default=1, help="Number of parallel processes (default=1)")
        return parser.parse_args()
    
def parse_pfamscan(file_content):
    """
    Parse pfamscan output content into a pandas DataFrame.
    Info: https://github.com/aziele/pfam_scan
    """
    lines = file_content.split('\n')
    data_lines = [
        line for line in lines if line and not line.startswith(('#', '#HMM', '#MATCH', '#PP', '#SEQ', '#CS'))
    ]
    # Define columns
    columns = [
        'seq_id', 'align_start', 'align_end', 'env_start', 'env_end',
        'hmm_acc', 'hmm_name', 'type', 'hmm_start', 'hmm_end',
        'hmm_length', 'bit_score', 'e_value', 'significance', 'clan'
    ]
    # Parse into DataFrame
    parsed_data = []
    for line in data_lines:
        values = line.split()
        if len(values) == len(columns):
            row = dict(zip(columns, values))
            for col in ['align_start', 'align_end', 'env_start', 'env_end', 'hmm_start', 'hmm_end', 'hmm_length']:
                row[col] = int(row[col])
            row['bit_score'] = float(row['bit_score'])
            row['e_value'] = float(row['e_value'])
            parsed_data.append(row)
    df = pd.DataFrame(parsed_data)
    df[['transcript_id', 'gene_id', 'otthumg', 'otthumt',
        'protein_name', 'gene_name', 'length']] = df['seq_id'].str.split('|', expand=True)
    return df


def hmmfetch_runner(hmm_database, hmm_id):
    """
    Fetches a specific HMM profile from an HMM database using the `hmmfetch` command.
    This function creates a temporary file to store the fetched HMM profile to a temporary file.
    """
    temp_file = tempfile.NamedTemporaryFile(delete=False, suffix=".hmm")
    output_file = temp_file.name
    temp_file.close()
    command = ["micromamba", "run", "hmmfetch", "-o", output_file, hmm_database, hmm_id]
    subprocess.run(command, check=True)
    return output_file


def hmmemit_runner(domain_id, input_hmm_file, output_fasta_file):
    """
    Computes consensus for a HMM profile.
    """
    command = ["micromamba", "run","hmmemit", "-o", output_fasta_file, "-C", input_hmm_file]
    try:
        subprocess.run(command, check=True)
        print(f"hmmemit executed successfully for {domain_id}")
    except subprocess.CalledProcessError as e:
        print(f"Error executing hmmemit for domain {domain_id}: {e}")


def concat_hmm_functions(hmmalign_output, consensus_file):
    """
    Concatenate the output of hmmalign and the consensus sequence.
    """
    output_file = consensus_file + ".with.consensus.fa"
    with open(output_file, 'wb') as outfile:
        for file in [hmmalign_output, consensus_file]:
            with open(file, 'rb') as infile:
                shutil.copyfileobj(infile, outfile)
    return output_file


def hmmalign_runner(input_fasta, input_hmm, output_file):
    """
    Performs sequence alignment using file using the hmmalign command.
    """
    command = [
        "micromamba", "run","hmmalign", "-o", output_file, "--amino", '--outformat', 'Pfam', input_hmm, input_fasta
    ]
    try:
        subprocess.run(command, check=True)
        print(f"hmmalign executed successfully: {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error executing hmmalign: {e}")
    return output_file


def extract_domain_sequences(fasta_file, pfam_df, padding):
    """
    Extracts domain sequences from the FASTA file based on the Pfam DataFrame.
    """
    domain_sequences = {}
    if fasta_file.endswith('.gz'):
        with gzip.open(fasta_file, "rt") as handle:
            sequences = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    else:
        with open(fasta_file, "r") as handle:
            sequences = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    for _, row in pfam_df.iterrows():
        domain_id = row['domain_id']
        seq_id = row['seq_id']
        start = max(0, row['align_start'] - padding)
        end = min(len(sequences[seq_id].seq), row['align_end'] + padding)
        domain_sequences[f"{seq_id}_{domain_id}_{start + 1}_{end}"] = sequences[seq_id].seq[start:end]
    return domain_sequences


def save_domain_sequences(domain_sequences, output_file):
    """
    Save the extracted domain sequences to a FASTA file.
    """
    with open(output_file, 'w') as f:
        for domain_id, sequence in domain_sequences.items():
            f.write(f">{domain_id}\n{sequence}\n")

def process_domain(hmm_acc, pfam_df, fasta_path, padding, output_dir, hmmdatabase):
    """
    Process a single domain by fetching the HMM profile, extracting sequences, and aligning them.
    """
    domain_df = pfam_df[pfam_df['hmm_acc'] == hmm_acc]
    print(f"[PID {os.getpid()}] Processing domain: {hmm_acc}")

    try:
        hmm_file = hmmfetch_runner(hmmdatabase, hmm_acc)
        consensus_file = tempfile.NamedTemporaryFile(delete=False, mode="w").name
        hmmemit_runner(hmm_acc, hmm_file, consensus_file)

        domain_sequences = extract_domain_sequences(fasta_path, domain_df, padding)
        temp_fasta = tempfile.NamedTemporaryFile(delete=False, mode="w").name
        save_domain_sequences(domain_sequences, temp_fasta)

        concatenated_file = concat_hmm_functions(temp_fasta, consensus_file)
        output_alignment = os.path.join(output_dir, f"{hmm_acc}.hmmalign.aln")
        hmmalign_runner(concatenated_file, hmm_file, output_alignment)

    finally:
        for f in [temp_fasta, consensus_file, hmm_file]:
            if os.path.exists(f):
                os.unlink(f)


def main():
        
    args = parse_arguments()

    with open(args.pfamscan, 'r') as f:
        pfam_content = f.read()
    pfam_df = parse_pfamscan(pfam_content)
    pfam_df['domain_id'] = pfam_df['hmm_acc'].str.split('.', expand=True)[0]

    os.makedirs(args.output_dir, exist_ok=True)

    hmm_accs = pfam_df['hmm_acc'].unique()

    with ProcessPoolExecutor(max_workers=args.n_cores) as executor:
        futures = [
            executor.submit(process_domain, hmm_acc, pfam_df, args.fasta, args.padding, args.output_dir, args.hmmdatabase)
            for hmm_acc in hmm_accs
        ]
        for future in futures:
            future.result()  # Ensure all tasks complete or raise exceptions

if __name__ == "__main__":
    main()
