"""Quality control for the MetaDome data release.

All checks read only released files. None needs the database, the prebuild
output or the application, so anyone holding the Zenodo deposit can run them.

A. Aggregation statistics. The tolerance track carries one row per genomic
   codon, but several transcripts can cover a codon and report different
   scores, in which case the median is taken. This recomputes the statistics
   from the final dataset and verifies the shipped track against them.

B. Regression against MetaDome v1.0.1 (2022, Zenodo record 6625251). The v1
   table covers GRCh37 and uses the same tolerance definition, so the codon
   sets and the scores should largely coincide.

C. Cross-build comparison, keyed on UniProt accession and residue position,
   which carries no genomic coordinate and so aligns the two assemblies.
"""

import argparse
import csv
import gzip
import os
import sys
from statistics import median

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))

from scripts.data_release.generate_derived_tracks import (
    run_sort, SORT_CODON, ZENODO_CURRENT, ZENODO_PREVIOUS)

SORT_RESIDUE = ["-k1,1", "-k2,2n"]

# The 2022 table is comma-separated despite its .tsv.gz name. Its columns are
# the final dataset's, without refseq_ids and without the protein block.
V1_DELIMITER = ","

DNDS_TOLERANCE = 1e-9


def open_table(path, delimiter="\t"):
    handle = gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") \
        else open(path, newline="", encoding="utf-8")
    return handle, csv.DictReader(handle, delimiter=delimiter)


def normalise(path, delimiter, destination):
    """Reduce a dataset to one row per (codon, transcript).

    A codon repeats once per domain placement with identical tolerance each
    time. Those repeats are consecutive, so comparing against the previous key
    collapses them without holding millions of keys in memory.
    """
    handle, reader = open_table(path, delimiter)
    previous, kept = None, 0
    with handle, open(destination, "w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        for row in reader:
            key = (row["chrom"], row["pos_start"], row["pos_stop"],
                   row["strand"], row["gencode_transcription_id"])
            if key == previous:
                continue
            previous = key
            writer.writerow([row["chrom"], row["pos_start"], row["pos_stop"],
                             row["strand"], row["sw_dn_ds"], row["sw_coverage"],
                             row["sw_size"], row["symbol"]])
            kept += 1
    return kept


def codon_view(path, delimiter, workdir, label):
    """Sorted one-row-per-codon view, plus the disagreement statistics."""
    unsorted = os.path.join(workdir, label + ".unsorted")
    ordered = os.path.join(workdir, label + ".sorted")
    rows = normalise(path, delimiter, unsorted)
    run_sort(unsorted, ordered, SORT_CODON)
    os.remove(unsorted)

    stats = dict(rows=rows, groups=0, gene_overlap=0,
                 size_mismatch=0, coverage_mismatch=0, dnds_mismatch=0)
    collapsed = os.path.join(workdir, label + ".codons")

    def flush(group, writer):
        if not group:
            return
        if not any(g[4] for g in group):
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
        writer.writerow([group[0][0], group[0][1], group[0][2], group[0][3],
                         repr(median(dnds)), repr(median(coverage))])

    with open(ordered, newline="", encoding="utf-8") as src, \
            open(collapsed, "w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        current, group = None, []
        for row in csv.reader(src, delimiter="\t"):
            key = (row[0], row[1], row[2], row[3])
            if key != current:
                flush(group, writer)
                current, group = key, []
            group.append(row)
        flush(group, writer)
    os.remove(ordered)
    return collapsed, stats


def residue_view(path, delimiter, workdir, label):
    """Sorted one-row-per-residue view, keyed on UniProt accession and position.

    That key carries no genomic coordinate, so it aligns the two assemblies
    directly. A residue repeats per transcript and per domain placement; the
    median is taken where transcripts disagree, as in the codon view.
    """
    unsorted = os.path.join(workdir, label + ".residues.unsorted")
    ordered = os.path.join(workdir, label + ".residues.sorted")
    handle, reader = open_table(path, delimiter)
    previous = None
    with handle, open(unsorted, "w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        for row in reader:
            if not row["protein_ac"] or row["sw_dn_ds"] == "":
                continue
            key = (row["protein_ac"], row["protein_pos"],
                   row["gencode_transcription_id"])
            if key == previous:
                continue
            previous = key
            writer.writerow([row["protein_ac"], row["protein_pos"], row["sw_dn_ds"]])
    run_sort(unsorted, ordered, SORT_RESIDUE)
    os.remove(unsorted)

    collapsed = os.path.join(workdir, label + ".residues")
    count = 0
    with open(ordered, newline="", encoding="utf-8") as src, \
            open(collapsed, "w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        current, values = None, []
        for row in csv.reader(src, delimiter="\t"):
            key = (row[0], row[1])
            if key != current:
                if values:
                    writer.writerow([current[0], current[1], repr(median(values))])
                    count += 1
                current, values = key, []
            values.append(float(row[2]))
        if values:
            writer.writerow([current[0], current[1], repr(median(values))])
            count += 1
    os.remove(ordered)
    return collapsed, count


def read_track(path):
    """Yield (key, dn_ds, coverage) from a released tolerance track.

    The track is 0-based half-open BED file; the final dataset is
    1-based, so the start is shifted back here to compare the two.
    """
    handle = gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") \
        else open(path, encoding="utf-8")
    with handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, start, stop, strand = fields[0], fields[1], fields[2], fields[5]
            yield (chrom, int(start) + 1, int(stop), strand), \
                float(fields[6]), float(fields[7])


def read_codons(path):
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            yield (row[0], int(row[1]), int(row[2]), row[3]), \
                float(row[4]), float(row[5])


def read_residues(path):
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            yield (row[0], int(row[1])), float(row[2]), 0.0


def compare(left, right, tolerance):
    """Merge-join two sorted streams. Both are in LC_ALL=C order."""
    result = dict(left_only=0, right_only=0, shared=0, agree=0)
    a = next(left, None)
    b = next(right, None)
    while a is not None and b is not None:
        if a[0] == b[0]:
            result["shared"] += 1
            if abs(a[1] - b[1]) <= tolerance:
                result["agree"] += 1
            a, b = next(left, None), next(right, None)
        elif a[0] < b[0]:
            result["left_only"] += 1
            a = next(left, None)
        else:
            result["right_only"] += 1
            b = next(right, None)
    while a is not None:
        result["left_only"] += 1
        a = next(left, None)
    while b is not None:
        result["right_only"] += 1
        b = next(right, None)
    return result


def percent(part, whole):
    return "{:.2f}%".format(100.0 * part / whole) if whole else "n/a"


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--final-dataset", required=True)
    parser.add_argument("--tolerance-track", required=True)
    parser.add_argument("--genome-build", required=True)
    parser.add_argument("--v1-table", help="2022 dn/ds table; enables check B")
    parser.add_argument("--cross-build-dataset", required=True,
                        help="final dataset for the other assembly; enables check C")
    parser.add_argument("--cross-build-name")
    parser.add_argument("--workdir", default=".")
    parser.add_argument("--report", help="Markdown output; stdout when omitted")
    parser.add_argument("--previous-track",
                        help="tolerance track from an earlier release; enables check D")
    args = parser.parse_args()

    lines = ["# MetaDome data release quality control",
             "",
             "Verification of the MetaDome {} data release. Each check below reads only "
             "the files listed in the table and needs no database, prebuild output or "
             "running application, so it can be reproduced from the Zenodo record "
             "alone.".format(args.genome_build),
             "",
             "Release record: <{}>".format(ZENODO_CURRENT),
             "Previous release: <{}>".format(ZENODO_PREVIOUS),
             "",
             "| input | file |",
             "| --- | --- |",
             "| final dataset | `{}` |".format(os.path.basename(args.final_dataset)),
             "| tolerance track | `{}` |".format(os.path.basename(args.tolerance_track))]
    if args.v1_table:
        lines.append("| v1.0.1 table | `{}` |".format(os.path.basename(args.v1_table)))
    if args.cross_build_dataset:
        lines.append("| cross-build dataset | `{}` |".format(os.path.basename(args.cross_build_dataset)))
    if args.previous_track:
        against = compare(read_track(args.previous_track),
                          read_track(args.tolerance_track), DNDS_TOLERANCE)
        union = against["shared"] + against["left_only"] + against["right_only"]
        lines.append("## Comparison with an earlier release")
        lines.append("")
        lines.append("| statistic | value |")
        lines.append("| --- | --- |")
        lines.append("| codons in `{}` | {} |".format(
            os.path.basename(args.previous_track),
            against["shared"] + against["left_only"]))
        lines.append("| codons in this release | {} |".format(
            against["shared"] + against["right_only"]))
        lines.append("| shared | {} ({} of all codons in either) |".format(
            against["shared"], percent(against["shared"], union)))
        lines.append("| earlier release only | {} |".format(against["left_only"]))
        lines.append("| this release only | {} |".format(against["right_only"]))
        lines.append("| shared codons with dn/ds equal to within {} | {} ({}) |".format(
            DNDS_TOLERANCE, against["agree"], percent(against["agree"], against["shared"])))
        lines.append("")
    lines.append("")

    codons, stats = codon_view(args.final_dataset, "\t", args.workdir, "v2")

    lines.append("## {} aggregation statistics".format(args.genome_build))
    lines.append("")
    lines.append("| statistic | value | share of codons |")
    lines.append("| --- | --- | --- |")
    lines.append("| codon and transcript combinations | {} | |".format(stats["rows"]))
    lines.append("| unique genomic codons | {} | |".format(stats["groups"]))
    for key, name in (("gene_overlap", "covered by more than one gene"),
                      ("size_mismatch", "transcripts disagree on window size"),
                      ("coverage_mismatch", "transcripts disagree on coverage"),
                      ("dnds_mismatch", "transcripts disagree on dn/ds")):
        lines.append("| {} | {} | {} |".format(
            name, stats[key], percent(stats[key], stats["groups"])))
    lines.append("")

    shipped = compare(read_codons(codons), read_track(args.tolerance_track), 0.0)
    consistent = (shipped["left_only"] == 0 and shipped["right_only"] == 0
                  and shipped["agree"] == shipped["shared"])
    lines.append("The tolerance track reproduces from the final dataset: "
                 "**{}** ({} codons compared, {} identical).".format("yes" if consistent else "no", shipped["shared"], shipped["agree"]))
    lines.append("")

    if args.v1_table:
        v1_codons, v1_stats = codon_view(args.v1_table, V1_DELIMITER, args.workdir, "v1")
        against = compare(read_codons(v1_codons), read_codons(codons), DNDS_TOLERANCE)
        union = against["shared"] + against["left_only"] + against["right_only"]
        lines.append("## Comparison with MetaDome v1.0.1 (2022)")
        lines.append("")
        lines.append("| statistic | value |")
        lines.append("| --- | --- |")
        lines.append("| codons in v1.0.1 | {} |".format(v1_stats["groups"]))
        lines.append("| codons in this release | {} |".format(stats["groups"]))
        lines.append("| shared | {} ({} of all codons in either release) |".format(
            against["shared"], percent(against["shared"], union)))
        lines.append("| v1.0.1 only | {} |".format(against["left_only"]))
        lines.append("| this release only | {} |".format(against["right_only"]))
        lines.append("| shared codons with dn/ds equal to within {} | {} ({}) |".format(
            DNDS_TOLERANCE, against["agree"], percent(against["agree"], against["shared"])))
        lines.append("")
        os.remove(v1_codons)

    if args.cross_build_dataset:
        here, here_count = residue_view(args.final_dataset, "\t", args.workdir, "this")
        there, there_count = residue_view(args.cross_build_dataset, "\t", args.workdir, "other")
        across = compare(read_residues(here), read_residues(there), DNDS_TOLERANCE)
        union = across["shared"] + across["left_only"] + across["right_only"]
        lines.append("## {} against {}".format(args.genome_build, args.cross_build_name))
        lines.append("")
        lines.append("Residues are matched by UniProt accession and position, which is "
                     "independent of the genome assembly. The two builds annotate "
                     "different GENCODE and gnomAD releases, so the final row is "
                     "expected to be small: identical tolerance requires identical "
                     "variant content in the sliding window, which different gnomAD "
                     "releases rarely produce. It is not a measure of disagreement "
                     "between the builds.")
        lines.append("")
        lines.append("| statistic | value |")
        lines.append("| --- | --- |")
        lines.append("| residues in {} | {} |".format(args.genome_build, here_count))
        lines.append("| residues in {} | {} |".format(args.cross_build_name, there_count))
        lines.append("| shared | {} ({} of all residues in either build) |".format(
            across["shared"], percent(across["shared"], union)))
        lines.append("| {} only | {} |".format(args.genome_build, across["left_only"]))
        lines.append("| {} only | {} |".format(args.cross_build_name, across["right_only"]))
        lines.append("| shared residues with dn/ds equal to within {} | {} ({}) |".format(
            DNDS_TOLERANCE, across["agree"], percent(across["agree"], across["shared"])))
        lines.append("")
        os.remove(here)
        os.remove(there)

    os.remove(codons)

    report = "\n".join(lines) + "\n"
    if args.report:
        with open(args.report, "w", encoding="utf-8") as handle:
            handle.write(report)
        print("wrote {}".format(args.report))
    else:
        sys.stdout.write(report)


if __name__ == "__main__":
    main()
