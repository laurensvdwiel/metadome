import argparse
from multiprocess import Pool, Manager
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import os
from tqdm import tqdm
import csv
import gzip

global_swissprot = {}
global_gencode_prot = {}
global_gencode_trans = {}
global_transcript_lookup = {}


def parse_args():
    parser = argparse.ArgumentParser(description="Validate perfect matches between Gencode and SwissProt using BLAST results.")
    parser.add_argument('--blast', required=True, help='Path to the BLAST results file (CSV format).')
    parser.add_argument('--swissprot', required=True, help='Path to SwissProt FASTA file.')
    parser.add_argument('--gencode_prot', required=True, help='Path to Gencode protein FASTA file.')
    parser.add_argument('--gencode_trans', required=True, help='Path to Gencode protein coding transcript FASTA file.')
    parser.add_argument('--out_pass', default="perfect_matches.csv", help='Output file for perfect matches.')
    parser.add_argument('--out_fail', default="failed_matches.csv", help='Output file for failed matches.')
    parser.add_argument('--n_cores', type=int, default=os.cpu_count(), help='Number of cores to use for parallel processing.')
    return parser.parse_args()



def open_fasta(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def init_globals(swissprot, gencode_prot, gencode_trans, transcript_lookup):
    global global_swissprot, global_gencode_prot, global_gencode_trans, global_transcript_lookup
    global_swissprot = swissprot
    global_gencode_prot = gencode_prot
    global_gencode_trans = gencode_trans
    global_transcript_lookup = transcript_lookup


def interpret_blast(line):
    """
    Interpret a line from the BLAST results file and return a dictionary of relevant fields.
    """
    result_line = line.strip().split(',')
    if len(result_line) < 12:
        raise ValueError("Malformed BLAST line with fewer than 12 fields")
    return {
        'qseqid': result_line[0],
        'sseqid': result_line[1],
        'pident': float(result_line[2]),
        'length': int(result_line[3]),
        'mismatch': int(result_line[4]),
        'gapopen': int(result_line[5]),
        'qstart': int(result_line[6]),
        'qend': int(result_line[7]),
        'sstart': int(result_line[8]),
        'send': int(result_line[9]),
        'evalue': float(result_line[10]),
        'bitscore': float(result_line[11])
    }


def process_line(line_data):
    line_number, line = line_data
    try:
        blast_result = interpret_blast(line)
        qseqid = blast_result['qseqid']
        sseqid = blast_result['sseqid']
        swissprot_id = sseqid.split('|')[1] if '|' in sseqid else sseqid

        swissprot_seq_records = [r for k, r in global_swissprot.items() if str("|" + swissprot_id + "|") in k]
        gencode_prot_seq_record = global_gencode_prot.get(qseqid, None)

        results = []

        if not swissprot_seq_records or not gencode_prot_seq_record:
            return [('fail', (qseqid, sseqid, "Missing protein sequence"))]

        for swissprot_seq_record in swissprot_seq_records:
            if str(swissprot_seq_record.seq) != str(gencode_prot_seq_record.seq):
                results.append(('fail', (qseqid, sseqid, "Sequence mismatch")))
                continue
            transcript_id = qseqid.split('|')[0] if '|' in qseqid else None

            if not transcript_id:
                results.append(('fail', (qseqid, sseqid, "Missing transcript ID in qseqid")))
                continue

            matching_key = next((tid for tid in global_transcript_lookup.keys() if tid.split('|')[0] == transcript_id), None)

            if not matching_key:
                results.append(('fail', (qseqid, sseqid, f"Transcript ID {transcript_id} matching key {matching_key} not found")))
                continue

            transcript_record = global_transcript_lookup[matching_key]

            try:
                cds_field = next((field for field in matching_key.split('|') if field.startswith("CDS:")), None)

                if not cds_field:
                    results.append(('fail', (qseqid, sseqid, f"Missing CDS information in transcript ID: {matching_key} transcript_record {transcript_record}")))
                    continue

                cds_start, cds_end = map(int, cds_field.replace("CDS:", "").split("-"))

            except (IndexError, ValueError) as e:
                results.append(('fail', (qseqid, sseqid, f"Failed to extract CDS start/end positions: {e}")))
                continue

            transcript_seq = transcript_record.seq[cds_start - 1:cds_end]

            if len(transcript_seq) % 3 != 0:
                results.append(('fail', (qseqid, sseqid, "Transcript length not a multiple of 3")))
                continue

            # if not transcript_seq.startswith(["ATG", ""]):
                # results.append(('fail', (qseqid, sseqid, "Transcript missing start codon")))
                # continue

            if transcript_seq[-3:] not in ["TAA", "TAG", "TGA"]:
                results.append(('fail', (qseqid, sseqid, "Transcript missing stop codon")))
                continue

            translated = Seq(str(transcript_seq)).translate(to_stop=True)
            if str(translated) != str(gencode_prot_seq_record.seq):
                results.append(('fail', (qseqid, sseqid, "Translation mismatch")))
                continue

            results.append(('pass', (qseqid, sseqid, transcript_record.id, swissprot_seq_record.id)))

        return results

    except Exception as e:
        return [('fail', (line_data[1].strip().split(',')[0], '', f"Exception: {str(e)}"))]


def write_result(result, pass_writer, fail_writer, f_pass, f_fail):
    for status, row in result:
        if status == 'pass':
            pass_writer.writerow(row)
        else:
            fail_writer.writerow(row)
    f_pass.flush()
    f_fail.flush()



def main():
    args = parse_args()

    print("Reading SwissProt FASTA...")
    with open_fasta(args.swissprot) as handle:
        swissprot = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

    print("Reading Gencode protein FASTA...")
    with open_fasta(args.gencode_prot) as handle:
        gencode_prot = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

    print("Reading Gencode transcript FASTA...")
    with open_fasta(args.gencode_trans) as handle:
        gencode_trans = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

    transcript_lookup = {}
    for k, v in gencode_trans.items():
        tid = k
        transcript_lookup[tid] = v

    print(f"Reading BLAST results from {args.blast}...")
    with open(args.blast, 'r') as file:
        lines = list(enumerate(file, start=1))

    print(f"Starting parallel processing with {args.n_cores} cores...")

    with open(args.out_pass, 'w', newline='') as f_pass, open(args.out_fail, 'w', newline='') as f_fail:
        pass_writer = csv.writer(f_pass)
        fail_writer = csv.writer(f_fail)

        # Write headers
        pass_writer.writerow(["qseqid", "sseqid", "transcript_id", "swissprot_id"])
        fail_writer.writerow(["qseqid", "sseqid", "reason"])
        f_pass.flush()
        f_fail.flush()

        with Pool(processes=args.n_cores, initializer=init_globals,
                initargs=(swissprot, gencode_prot, gencode_trans, transcript_lookup)) as pool:
            for result in tqdm(pool.imap_unordered(process_line, lines), total=len(lines)):
                write_result(result, pass_writer, fail_writer, f_pass, f_fail)

        print(f"Perfect matches written to {args.out_pass}")
        print(f"Failed matches written to {args.out_fail}")
        print("Processing complete.")
    
    

if __name__ == "__main__":
    main()