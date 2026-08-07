"""Derive the release track files from the MetaDome final dataset.

Each track is built from the dataset produced by generate_final_dataset.py, so
they always describe the same data.

    tolerance-landscape   dn/ds over a 21-codon sliding window, one row per codon
    metadomain-clinvar    pathogenic ClinVar projected through Pfam alignments
    pfam-domain-coverage  every codon aligning to a Pfam domain consensus position

Genomic positions are 1-based inclusive. Where transcripts overlap a codon and
their tolerance scores disagree, the median is used.
"""

import argparse
import csv
import gzip
import os
import shutil
import subprocess

DATASET_COLUMNS = [
    "chrom", "pos_start", "pos_stop", "strand",
    "symbol", "gencode_transcription_id", "refseq_ids",
    "sw_dn_ds", "sw_coverage", "sw_size",
    "protein_ac", "protein_pos", "ref_aa", "ref_codon", "cdna_pos", "exon_numbers",
    "domain_id", "consensus_pos",
    "normal_variant_count", "normal_missense_variant_count",
    "pathogenic_variant_count", "pathogenic_missense_variant_count",
    "pathogenic_P_count", "pathogenic_LP_count",
    "pathogenic_missense_P_count", "pathogenic_missense_LP_count",
    "meta_domain_clinvar_P_records", "meta_domain_clinvar_LP_records",
]

# The pathogenic and likely-pathogenic colours used in the MetaDome web interface.
COLOUR_PATHOGENIC = "214,39,40"        # #d62728
COLOUR_LIKELY_PATHOGENIC = "255,127,14"  # #ff7f0e

# Embedded in every feature as a link back to MetaDome.
BASE_URL = "https://www.metadome.app/metadome"

# Zenodo record for this release, and the one preceding it. The 2022 record
# holds the data for the MetaDome paper, doi:10.1002/humu.23798.
ZENODO_CURRENT = "https://zenodo.org/records/19376150"
ZENODO_PREVIOUS = "https://zenodo.org/records/6625251"

def open_final_dataset(path):
    handle = gzip.open(path, "rt", newline="", encoding="utf-8") if path.endswith(".gz") \
        else open(path, newline="", encoding="utf-8")
    reader = csv.DictReader(handle, delimiter="\t")
    missing = [c for c in DATASET_COLUMNS if c not in (reader.fieldnames or [])]
    if missing:
        raise SystemExit("Final dataset file is missing columns: {}".format(", ".join(missing)))
    return handle, reader


def median(values):
    if not values:
        return ""
    values = sorted(values)
    middle = len(values) // 2
    if len(values) % 2:
        return values[middle]
    return (values[middle - 1] + values[middle]) / 2.0


def _has(tool):
    from shutil import which
    return which(tool) is not None


SORT_BED = ["-k1,1", "-k2,2n"]                          # what bedToBigBed requires
SORT_CODON = ["-k1,1", "-k2,2n", "-k3,3n", "-k4,4"]     # full grouping key

def run_sort(source, destination, keys=SORT_BED):
    """Sort a TSV. `keys` must cover the whole grouping key when the caller then
    streams adjacent rows into groups - sorting on only chrom and start leaves
    rows that share those two in unspecified relative order, which would split
    one codon's rows across several groups.

    LC_ALL=C is not optional: under a locale-aware collation the chrom ordering
    differs and bedToBigBed rejects the result. GNU sort is preferred when
    present because BSD sort has no --parallel.
    """
    tool = "gsort" if _has("gsort") else "sort"
    command = [tool, "-t", "\t"] + list(keys) + ["-S", "2G"]
    if tool == "gsort":
        command += ["--parallel", "4"]
    command += [source, "-o", destination]
    subprocess.check_call(command, env=dict(os.environ, LC_ALL="C"))


def header_block(prefix, genome_build, role, notes=()):
    """Commented provenance block for a release file.

    The release stem is embedded verbatim rather than reconstructed, so the
    header cannot drift from the filename.
    """
    lines = ["MetaDome data release - " + role,
             "assembly: " + genome_build,
             "release: " + prefix]
    lines.extend(notes)
    lines.append("coordinates: 1-based inclusive")
    lines.append("release record: " + ZENODO_CURRENT)
    lines.append("previous release: " + ZENODO_PREVIOUS)
    lines.append("https://www.metadome.app/metadome")
    return "".join("# " + line + "\n" for line in lines)


def prepend_header(path, header):
    """Add the header after sorting; sorting would scatter the comment lines."""
    staged = path + ".hdr"
    with open(staged, "w", encoding="utf-8") as out:
        out.write(header)
        with open(path, encoding="utf-8") as body:
            shutil.copyfileobj(body, out)
    os.replace(staged, path)

# --------------------------------------------------------------------------
# B: tolerance BED
# --------------------------------------------------------------------------

def build_tolerance_bed(final_dataset, destination, header):
    """Deduplicate to one row per genomic codon, median where transcripts disagree.

    The dataset is grouped by (chrom, pos_start, pos_stop, strand). Rows for the
    same codon are adjacent only after sorting, so this streams a pre-sorted
    intermediate rather than holding ~11M groups in memory.
    """
    intermediate = destination + ".unsorted"
    handle, reader = open_final_dataset(final_dataset)
    previous = None
    with handle, open(intermediate, "w", newline="", encoding="utf-8") as tmp:
        writer = csv.writer(tmp, delimiter="\t", lineterminator="\n")
        for row in reader:
            # A codon repeats once per domain placement, with identical
            # tolerance each time. Those repeats are written consecutively by
            # generate_final_dataset, so comparing with the previous key collapses
            # them in constant memory - a `seen` set would grow to ~24M entries.
            key = (row["chrom"], row["pos_start"], row["pos_stop"], row["strand"],
                   row["gencode_transcription_id"])
            if key == previous:
                continue
            previous = key
            writer.writerow([row["chrom"], row["pos_start"], row["pos_stop"], row["strand"],
                             row["sw_dn_ds"], row["sw_coverage"], row["sw_size"], row["symbol"]])

    sorted_path = destination + ".sorted"
    run_sort(intermediate, sorted_path, SORT_CODON)
    os.remove(intermediate)

    stats = dict(rows=0, groups=0, gene_overlap=0,
                 size_mismatch=0, coverage_mismatch=0, dnds_mismatch=0)

    def flush(group, out):
        if not group:
            return
        stats["groups"] += 1
        dnds = [float(g[4]) for g in group if g[4] != ""]
        coverage = [float(g[5]) for g in group if g[5] != ""]
        if len({g[6] for g in group}) > 1:
            stats["size_mismatch"] += 1
        if len(set(coverage)) > 1:
            stats["coverage_mismatch"] += 1
        if len(set(dnds)) > 1:
            stats["dnds_mismatch"] += 1
        if len({g[7] for g in group}) > 1:
            stats["gene_overlap"] += 1
        info = "sw_dn_ds:{},sw_coverage:{}".format(median(dnds), median(coverage))
        out.writerow([group[0][0], group[0][1], group[0][2], group[0][3], info])

    with open(sorted_path, newline="", encoding="utf-8") as src, \
            open(destination, "w", newline="", encoding="utf-8") as out:
        out.write(header)
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["#chrom", "pos_start", "pos_stop", "strand", "info"])
        current_key, group = None, []
        for row in csv.reader(src, delimiter="\t"):
            stats["rows"] += 1
            key = (row[0], row[1], row[2], row[3])
            if key != current_key:
                flush(group, writer)
                current_key, group = key, []
            group.append(row)
        flush(group, writer)
    os.remove(sorted_path)
    return stats


# --------------------------------------------------------------------------
# C and D: the bigBed inputs
# --------------------------------------------------------------------------

def metadome_url(genome_build, chrom, position):
    return "{}/position/{}/{}/{}/".format(BASE_URL, genome_build, chrom, position)


def build_domain_tracks(final_dataset, genome_build, clinvar_path, coverage_path, clinvar_header, coverage_header):
    """C keeps codons in a Pfam domain that carry at least one homologous P/LP.
    D keeps every codon in a Pfam domain. Both emit one row per placement."""
    handle, reader = open_final_dataset(final_dataset)
    clinvar_rows = coverage_rows = 0

    with handle, \
            open(clinvar_path, "w", newline="", encoding="utf-8") as clinvar_out, \
            open(coverage_path, "w", newline="", encoding="utf-8") as coverage_out:
        clinvar_writer = csv.writer(clinvar_out, delimiter="\t", lineterminator="\n")
        coverage_writer = csv.writer(coverage_out, delimiter="\t", lineterminator="\n")

        for row in reader:
            domain = row["domain_id"]
            if not domain:
                continue

            chrom, start, stop = row["chrom"], row["pos_start"], row["pos_stop"]
            strand, consensus = row["strand"], row["consensus_pos"]
            accession, protein_pos = row["protein_ac"], row["protein_pos"]
            name = "{}/{}:{}:{}".format(accession, protein_pos, domain, consensus)
            url = metadome_url(genome_build, chrom, start)

            coverage_writer.writerow([
                chrom, start, stop, name, 0, strand,
                accession, protein_pos, domain, consensus, url,
            ])
            coverage_rows += 1

            pathogenic = int(row["pathogenic_P_count"] or 0)
            likely = int(row["pathogenic_LP_count"] or 0)
            if pathogenic + likely == 0:
                continue

            clinvar_writer.writerow([
                chrom, start, stop, name,
                min(1000, (pathogenic + likely) * 10),
                strand, start, stop,
                COLOUR_PATHOGENIC if pathogenic >= likely else COLOUR_LIKELY_PATHOGENIC,
                accession, protein_pos, domain, consensus,
                pathogenic, likely,
                row["pathogenic_missense_P_count"] or 0,
                row["pathogenic_missense_LP_count"] or 0,
                row["meta_domain_clinvar_P_records"],
                row["meta_domain_clinvar_LP_records"],
                url,
            ])
            clinvar_rows += 1

    for path, header in ((clinvar_path, clinvar_header), (coverage_path, coverage_header)):
          sorted_path = path + ".sorted"
          run_sort(path, sorted_path)
          os.replace(sorted_path, path)
          prepend_header(path, header)

    return clinvar_rows, coverage_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--final-dataset", required=True, help="final dataset .tsv or .tsv.gz")
    parser.add_argument("--genome-build", required=True, help="e.g. GRCh37.p13")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--prefix", required=True,
                        help="Filename stem, e.g. MetaDome_v2.0_GRCh37.p13_GENCODE-v19_...")
    parser.add_argument("--skip-tolerance", action="store_true")
    parser.add_argument("--skip-domains", action="store_true")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    stem = os.path.join(args.out_dir, args.prefix)

    if not args.skip_tolerance:
        destination = stem + "_derived-track-tolerance-landscape-sw10.bed"
        header = header_block(
            args.prefix, args.genome_build,
            "tolerance landscape, dn/ds over a 21-codon sliding window",
            ["one row per genomic codon",
             "median taken across transcripts that disagree"])
        print("building tolerance BED -> {}".format(destination), flush=True)
        stats = build_tolerance_bed(args.final_dataset, destination, header)
        total = stats["groups"] or 1
        print("  rows in            : {}".format(stats["rows"]))
        print("  unique codons out  : {}".format(stats["groups"]))
        print("  gene overlap       : {} ({:.2f}%)".format(
            stats["gene_overlap"], 100.0 * stats["gene_overlap"] / total))
        print("  sw_size mismatch   : {} ({:.2f}%)".format(
            stats["size_mismatch"], 100.0 * stats["size_mismatch"] / total))
        print("  coverage mismatch  : {} ({:.2f}%)".format(
            stats["coverage_mismatch"], 100.0 * stats["coverage_mismatch"] / total))
        print("  dn/ds mismatch     : {} ({:.2f}%)".format(
            stats["dnds_mismatch"], 100.0 * stats["dnds_mismatch"] / total))
        print("  (median applied where transcripts disagree)")
        print("  this release     : {}".format(ZENODO_CURRENT))
        print("  previous release : {}".format(ZENODO_PREVIOUS))

    if not args.skip_domains:
        clinvar_path = stem + "_derived-track-metadomain-clinvar.bed"
        coverage_path = stem + "_derived-track-pfam-domain-coverage.bed"
        print("building domain tracks -> {} , {}".format(clinvar_path, coverage_path), flush=True)
        clinvar_header = header_block(
            args.prefix, args.genome_build, "meta-domain ClinVar track, bed9+11",
            ["codons in a Pfam domain with at least one homologous pathogenic",
             "or likely pathogenic ClinVar variant at an equivalent position",
             "counts are of distinct variants; accession lists can be shorter",
             "field definitions in metadomain_clinvar.as"])
        coverage_header = header_block(
            args.prefix, args.genome_build, "Pfam domain coverage track, bed6+5",
            ["every codon aligning to a Pfam domain consensus position",
             "field definitions in pfam_coverage.as"])
        clinvar_rows, coverage_rows = build_domain_tracks(
            args.final_dataset, args.genome_build, clinvar_path, coverage_path,
            clinvar_header, coverage_header)
        print("  metadomain ClinVar rows : {}".format(clinvar_rows))
        print("  Pfam coverage rows      : {}".format(coverage_rows))
        # bed9 + 11 extras: uniprot_ac, uniprot_pos, pfam_id, consensus_pos,
        # four homologue counts, two accession lists, metadome_url.
        print("\nConvert with bedToBigBed (the comment header must be stripped):")
        print("  grep -v '^#' {} \\".format(clinvar_path))
        print("    | bedToBigBed -type=bed9+11 -as=metadomain_clinvar.as -tab \\")
        print("        stdin <chrom.sizes> {}".format(clinvar_path.replace(".bed", ".bb")))
        print("  grep -v '^#' {} \\".format(coverage_path))
        print("    | bedToBigBed -type=bed6+5 -as=pfam_coverage.as -tab \\")
        print("        stdin <chrom.sizes> {}".format(coverage_path.replace(".bed", ".bb")))


if __name__ == "__main__":
    main()
