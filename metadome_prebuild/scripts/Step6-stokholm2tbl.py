import os
import argparse
import pandas as pd
from Bio import AlignIO
from multiprocessing import Pool, cpu_count

def parse_arguments():
    parser = argparse.ArgumentParser(description="Compute metapositions from a folder of Stockholm alignments files")
    parser.add_argument("--input_folder", required=True, help="Path to folder containing Stockholm files.")
    parser.add_argument("--output", required=True, help="Path to output TSV file contain metapositions.")
    parser.add_argument("--n_cores", type=int, default=cpu_count(), help=f"Number of CPU cores to use (default = all available: {cpu_count()})")
    return parser.parse_args()

def map_alignment_positions(stockholm_file):
    """
    Map positions between sequences in a Stockholm alignment, using 'consensus' as reference.
    """
    try:
        alignment = AlignIO.read(stockholm_file, "stockholm")

        ref_id = [seq.id for seq in alignment if "consensus" in seq.id]
        if not ref_id:
            raise ValueError(f"No 'consensus' sequence found in {stockholm_file}")
        ref_seq = str(*[seq.seq for seq in alignment if "consensus" in seq.id])
        ref_index = [seq.id for seq in alignment].index(*ref_id)

        all_mappings = []

        for idx, record in enumerate(alignment):
            if idx == ref_index:
                continue

            seq_id = record.id
            seq = str(record.seq)

            try:
                transcript_id, domain_id, start, end = seq_id.split("_")
            except ValueError:
                raise ValueError(f"Invalid sequence ID format in {stockholm_file}: {seq_id}")

            ref_pos = 0
            seq_pos = int(start) - 1

            ref_positions = []
            seq_positions = []
            aa_ref = []
            aa_seq = []

            for i in range(len(ref_seq)):
                ref_char = ref_seq[i]
                seq_char = seq[i]

                if ref_char.upper().isalpha():
                    ref_pos += 1
                if seq_char.upper().isalpha():
                    seq_pos += 1

                ref_positions.append(ref_pos if ref_char.upper().isalpha() else None)
                seq_positions.append(seq_pos if seq_char.upper().isalpha() else None)
                aa_ref.append(ref_char)
                aa_seq.append(seq_char)

            df = pd.DataFrame({
                'transcript_id': transcript_id,
                'domain_id': domain_id,
                'range_id': f"{start}_{end}",
                'seq_position': seq_positions,
                'sequence_aa': aa_seq,
                'position_in_domain_consensus': ref_positions,
                'aa_in_domain_consensus': aa_ref
            })

            df = df.dropna(subset=["seq_position", "position_in_domain_consensus"])
            all_mappings.append(df)

        return pd.concat(all_mappings, ignore_index=True)

    except Exception as e:
        print(f"Error processing {stockholm_file}: {e}")
        return None


def process_stockholm_folder(input_folder, output_path, n_cores):
    """
    Process all Stockholm files in the input folder in parallel.
    """
    files = [
        os.path.join(input_folder, f)
        for f in os.listdir(input_folder)
        if os.path.isfile(os.path.join(input_folder, f))
    ]

    print(f"Found {len(files)} files. Using {n_cores} cores for processing...")

    with Pool(processes=n_cores) as pool:
        results = pool.map(map_alignment_positions, files)

    valid_results = [df for df in results if df is not None]

    if valid_results:
        final_df = pd.concat(valid_results, ignore_index=True)
        final_df.to_csv(output_path, sep="\t", index=False)
        print(f"\n All results written to: {output_path}")
    else:
        print("No valid alignments processed.")


def main():
    args = parse_arguments()
    print("Processing Stokcholm files")
    process_stockholm_folder(args.input_folder, args.output, args.n_cores)
    


if __name__ == "__main__":
    main()
