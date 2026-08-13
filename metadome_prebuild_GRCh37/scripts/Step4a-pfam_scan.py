import sys
import subprocess
import os
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Runs pfam_scan.pl on a target fasta file.")
    parser.add_argument("--input_pfam", required=True, help="Path to folder containing PFAM HMM.")
    parser.add_argument("--n_cores", type=int, help=f"Number of CPU cores to use.")
    parser.add_argument("--fasta_to_annotate", required=True, help="Path to the target protein FASTA file.")
    parser.add_argument("--manifest_path", required=True, help="Path to the toml.")
    parser.add_argument("--output_file", required=True, help="Path where to save pfamscan results.")

    return parser.parse_args()

def run_command(command, check=True):
    """Run a shell command and optionally check for errors."""
    print(f"Running: {' '.join(command)}")
    subprocess.run(command, check=check)

def main():
    args = parse_arguments()

    print("Provided variables are:")
    print(f"INPUT_PFAM: {args.input_pfam}")
    print(f"CPUs: {args.n_cores}")
    print(f"FASTA_TO_ANNOTATE: {args.fasta_to_annotate}")

    pfam_base_name = os.path.basename(args.input_pfam)

    pmf = args.manifest_path
    
    # # Preprocess Pfam HMMs
    # Needs directory cleanup otherwise it fails in case or resumed runs
    allowed_files = {"Pfam-A.hmm", "Pfam-A.hmm.dat"}
    for filename in os.listdir(args.input_pfam):
        if filename not in allowed_files:
            file_path = os.path.join(args.input_pfam, filename)
            if os.path.isfile(file_path):
                print(f"Deleting file: {file_path}")
                os.remove(file_path)
                
    print(f"Running hmmpress on {args.input_pfam}...")
    run_command(["pixi", "run", "--manifest-path", pmf , "hmmpress", os.path.join(args.input_pfam, "Pfam-A.hmm")])

    # Run Pfam Scan
    output_file_pfamscan = args.output_file
    print(f"Running pfam_scan.pl on {args.fasta_to_annotate}...")

    run_command([
        "pixi", "run","--manifest-path", pmf , "pfam_scan.pl",
        "-outfile", output_file_pfamscan,
        "-cpu", str(args.n_cores),
        "-fasta", args.fasta_to_annotate,
        "-dir", args.input_pfam
    ])

    print(f"Pfam Scan completed. Output saved to {output_file_pfamscan}.")

if __name__ == "__main__":
    main()
