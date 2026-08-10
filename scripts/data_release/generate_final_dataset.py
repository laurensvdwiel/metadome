"""Build the MetaDome final dataset for one genome build.

Reads prebuild output from disk. No database and no running application are
required.

  metadome_visualization/<build>/<transcript>/metadome_visualization.json
      codon position, strand, protein coordinates, tolerance, per-domain counts
  metadomains/<build>/<PFxxxxx>/metadomain_snv_annotation_clinvar
      ClinVar accessions, keyed by domain consensus position

The derived track files are built from this output by generate_derived_tracks.py.

Output is one row per transcript, codon fragment and domain placement. Codons
outside any Pfam domain are included with the domain columns empty, and a codon
split across an exon boundary yields one row per fragment.

Genomic positions are 1-based inclusive. Protein positions are 1-based from the
initiator methionine. The two sources are joined on (domain_id, consensus_pos).

Homologous counts and accessions describe the evolutionarily equivalent
positions in other genes. ClinVar records on the annotated codon itself are
reported separately in clinvar_P_at_position and clinvar_LP_at_position.
--check compares the counts against the variant counts recorded during the
prebuild and stops if they disagree.
"""

import argparse
import csv
import glob
import gzip
import multiprocessing
import os
import re
import shutil
import sys
import time
from multiprocessing import Pool

# Workers share the ClinVar lookup copy-on-write. The spawn start method, the
# default on macOS, would pickle a full copy into each worker instead.
try:
    multiprocessing.set_start_method("fork")
except RuntimeError:
    pass

try:
    import orjson

    def _load(fh):
        return orjson.loads(fh.read())

    _READ_MODE = "rb"
    _JSON = "orjson"
except ImportError:
    import json

    def _load(fh):
        return json.load(fh)

    _READ_MODE = "r"
    _JSON = "json (pip install orjson to speed up parsing)"


COLUMNS = [
    "chrom", "pos_start", "pos_stop", "strand",
    "symbol", "gencode_transcription_id", "refseq_ids",
    "sw_dn_ds", "sw_coverage", "sw_size",
    "protein_ac", "protein_pos", "ref_aa", "ref_codon", "cdna_pos", "exon_numbers",
    "domain_id", "consensus_pos",
    "normal_variant_count", "normal_missense_variant_count",
    "pathogenic_variant_count", "pathogenic_missense_variant_count",
    "pathogenic_P_count", "pathogenic_LP_count",
    "pathogenic_missense_P_count", "pathogenic_missense_LP_count",
    "clinvar_P_at_position", "clinvar_LP_at_position",
    "meta_domain_clinvar_P_records", "meta_domain_clinvar_LP_records",
    "meta_domain_clinvar_P_missense_records", "meta_domain_clinvar_LP_missense_records",
]

BLANK_PLACEMENT = [""] * 16  # domain_id .. LP records, for codons outside any domain

VISUALIZATION_FILE = "metadome_visualization.json"
CLINVAR_FILE = "metadomain_snv_annotation_clinvar"

# Set per worker by _init_worker; read-only thereafter.
_CLINVAR = {}
_STRICT = False


class ConsistencyError(RuntimeError):
    """Accession lists disagreed with the prebuilt counts."""


# --------------------------------------------------------------------------
# configuration
# --------------------------------------------------------------------------

def repo_root():
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def read_env(root):
    """Parse the repo-root .env - the same file docker compose substitutes
    from - so one invocation resolves correctly on the host or in a container."""
    values = {}
    path = os.path.join(root, ".env")
    if not os.path.isfile(path):
        return values
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, _, value = line.partition("=")
            value = value.strip()
            if value[:1] in ('"', "'"):
                closing = value.find(value[0], 1)
                value = value[1:closing] if closing > 0 else value[1:]
            else:
                # An unquoted value ends at the first whitespace-preceded '#',
                # which is how docker compose reads the same file.
                value = re.split(r"\s+#", value, maxsplit=1)[0].strip()
            values[key.strip()] = value
    return values


def resolve_data_dir(explicit):
    return explicit or read_env(repo_root()).get("METADOME_DIR") or "/usr/data"


# --------------------------------------------------------------------------
# ClinVar lookup
# --------------------------------------------------------------------------

def load_clinvar(metadomain_dir):
    """{(domain_id, consensus_pos): [(clinvar_id, is_pathogenic, is_missense, representation)]}

    The lookup is held in memory; the ClinVar annotations total roughly 80 MB
    per genome build.
    """
    lookup = {}
    rows = 0
    for path in sorted(glob.glob(os.path.join(metadomain_dir, "*", CLINVAR_FILE))):
        domain_id = os.path.basename(os.path.dirname(path))
        if os.path.getsize(path) < 8:  # a domain with no ClinVar holds a bare '""'
            continue
        with open(path, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            if not reader.fieldnames or "consensus_pos" not in reader.fieldnames:
                continue
            for row in reader:
                identifier = (row.get("clinvar_ID") or "").strip()
                if not identifier:
                    continue
                try:
                    consensus = int(row["consensus_pos"])
                except (KeyError, TypeError, ValueError):
                    continue
                lookup.setdefault((domain_id, consensus), []).append(
                    (identifier,
                     row.get("clinvar_clinsig") == "Pathogenic",
                     row.get("variant_type") == "missense",
                     row.get("unique_snv_str_representation", ""))
                )
                rows += 1
    return lookup, rows


def homologues(domain_id, consensus, own_codon):
    """Counts and ClinVar accessions at a single consensus position, in other genes.

    own_codon is the annotated codon's Codon.unique_str_representation, which
    identifies it by chromosome, region tuples and strand. A record whose
    representation contains it sits on that codon and is not a homologue of it.

    Counts are of distinct single nucleotide variants. One ClinVar record can
    appear as more than one variant when transcripts place it in different
    reading frames, so the accession lists can be shorter than the counts.
    """
    pathogenic, likely = {}, {}
    pathogenic_ids, likely_ids = {}, {}
    pathogenic_missense_ids, likely_missense_ids = {}, {}
    at_position_p, at_position_lp = set(), set()
    for identifier, is_pathogenic, is_missense, representation in _CLINVAR.get((domain_id, consensus), ()):
        if own_codon in representation:
            (at_position_p if is_pathogenic else at_position_lp).add(representation)
            continue
        if is_pathogenic:
            pathogenic[representation] = is_missense
            pathogenic_ids[identifier] = None
            if is_missense:
                pathogenic_missense_ids[identifier] = None
        else:
            likely[representation] = is_missense
            likely_ids[identifier] = None
            if is_missense:
                likely_missense_ids[identifier] = None
    return {
        "p": len(pathogenic),
        "lp": len(likely),
        "p_missense": sum(1 for missense in pathogenic.values() if missense),
        "lp_missense": sum(1 for missense in likely.values() if missense),
        "p_missense_ids": list(pathogenic_missense_ids),
        "lp_missense_ids": list(likely_missense_ids),
        "at_p": len(at_position_p),
        "at_lp": len(at_position_lp),
        "p_ids": list(pathogenic_ids),
        "lp_ids": list(likely_ids),
    }

# --------------------------------------------------------------------------
# phase 1: flatten one transcript
# --------------------------------------------------------------------------

def parse_positions(chr_positions):
    """'g.69091-69093' -> [(69091, 69093)].

    A codon split across an exon boundary arrives as
    'g.128407650, g.128408290-128408291' and yields one range per fragment.
    """
    ranges = []
    for token in str(chr_positions).split(","):
        token = token.strip()
        if token.startswith("g."):
            token = token[2:]
        if not token:
            continue
        bounds = token.split("-")
        try:
            start = int(bounds[0])
            stop = int(bounds[1]) if len(bounds) > 1 else start
        except (ValueError, IndexError):
            continue
        ranges.append((start, stop))
    return ranges


def flatten(path):
    """Return (rows, mismatch_count) for one transcript."""
    with open(path, _READ_MODE) as fh:
        data = _load(fh)

    transcript = data.get("transcript_id", "")
    symbol = data.get("gene_name", "")
    protein_ac = data.get("protein_ac", "")
    refseq = data.get("refseq_ids", "")

    rows, mismatches = [], 0

    for entry in data.get("positional_annotation", ()):
        ranges = parse_positions(entry.get("chr_positions", ""))
        if not ranges:
            continue

        # Mirrors Codon.unique_str_representation: chromosome, region tuples and
        # strand. Positions alone are not enough, since PAR1 gives chrX and chrY
        # identical coordinates.
        own_codon = "{}:[{}]::(Strand.{})".format(
            entry.get("chr", ""),
            ", ".join("({}, {})".format(start, stop) for start, stop in ranges),
            "plus" if entry.get("strand") == "+" else "minus")

        shared = [
            entry.get("chr", ""), None, None, entry.get("strand", ""),
            symbol, transcript, refseq,
            entry.get("sw_dn_ds", ""), entry.get("sw_coverage", ""), entry.get("sw_size", ""),
            protein_ac, entry.get("protein_pos", ""), entry.get("ref_aa", ""),
            entry.get("ref_codon", ""), entry.get("cdna_pos", ""), entry.get("exon_numbers", ""),
        ]

        placements = []
        for domain_id, domain in (entry.get("domains") or {}).items():
            # A null entry means the codon lies within the domain span but
            # aligns to no consensus column, so there is nothing to join on.
            if not domain:
                continue
            consensus_positions = domain.get("consensus_pos") or []
            if not consensus_positions:
                continue

            per_class = domain.get("pathogenic_variant_count_per_clinsig") or {}
            per_class_missense = domain.get("pathogenic_missense_variant_count_per_clinsig") or {}
            expected = {
                "p": per_class.get("Pathogenic", 0),
                "lp": per_class.get("Likely_pathogenic", 0),
                "p_missense": per_class_missense.get("Pathogenic", 0),
                "lp_missense": per_class_missense.get("Likely_pathogenic", 0),
            }

            # Records are attributed to the consensus position they were found
            # at; a residue aligning to several columns keeps them separate.
            per_consensus = [(consensus, homologues(domain_id, consensus, own_codon))
                             for consensus in consensus_positions]
            totals = {key: sum(found[key] for _, found in per_consensus)
                      for key in expected}

            if totals != expected:
                mismatches += 1
                if _STRICT:
                    raise ConsistencyError(
                        "{} {} consensus={}: summed {} but prebuilt counts {}".format(
                            transcript, domain_id, consensus_positions, totals, expected)
                    )

            domain_counts = [
                domain.get("normal_variant_count", 0),
                domain.get("normal_missense_variant_count", 0),
                domain.get("pathogenic_variant_count", 0),
                domain.get("pathogenic_missense_variant_count", 0),
            ]
            for consensus, found in per_consensus:
                placements.append(
                    [domain_id, consensus] + domain_counts + [
                        found["p"], found["lp"], found["p_missense"], found["lp_missense"],
                        found["at_p"], found["at_lp"],
                        ",".join(found["p_ids"]), ",".join(found["lp_ids"]),
                        ",".join(found["p_missense_ids"]), ",".join(found["lp_missense_ids"]),
                    ])

        if not placements:
            placements = [BLANK_PLACEMENT]

        for start, stop in ranges:
            for placement in placements:
                row = list(shared)
                row[1], row[2] = start, stop
                rows.append(row + placement)

    return rows, mismatches


def _init_worker(clinvar, strict):
    global _CLINVAR, _STRICT
    _CLINVAR = clinvar
    _STRICT = strict


def _write_shard(job):
    index, paths, shard_path = job
    written = mismatches = 0
    with gzip.open(shard_path, "wt", newline="", encoding="utf-8", compresslevel=6) as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        for path in paths:
            try:
                rows, bad = flatten(path)
            except ConsistencyError:
                raise  # --check must abort the run, not be logged and skipped
            except Exception as exc:  # an unreadable transcript is logged and skipped
                print("WARN {}: {}".format(path, exc), file=sys.stderr)
                continue
            writer.writerows(rows)
            written += len(rows)
            mismatches += bad
    return index, written, mismatches


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------

def chunk(items, count):
    """Round-robin so every shard mixes large and small transcripts."""
    buckets = [[] for _ in range(count)]
    for position, item in enumerate(items):
        buckets[position % count].append(item)
    return [bucket for bucket in buckets if bucket]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--genome-build", required=True, help="e.g. GRCh37.p13 or GRCh38.p14")
    parser.add_argument("--data-dir", default=None,
                        help="Overrides METADOME_DIR from the repo-root .env; falls back to /usr/data.")
    parser.add_argument("--out", required=True, help="Output .tsv.gz path")
    parser.add_argument("--workers", type=int, default=4,
                        help="Keep below your core count while a prebuild is running (default 4).")
    parser.add_argument("--limit", type=int, default=None,
                        help="Only process the first N transcripts. Use for a timed pilot.")
    parser.add_argument("--check", action="store_true",
                        help="Abort on the first row whose accession lists disagree with the "
                             "prebuilt counts. Recommended for the pilot.")
    parser.add_argument("--report", action="store_true",
                        help="Report what would be processed, then exit without writing.")
    args = parser.parse_args()

    data_dir = resolve_data_dir(args.data_dir)
    visualization_dir = os.path.join(data_dir, "metadome_visualization", args.genome_build)
    # metadomains/ is keyed by assembly without the patch suffix: GRCh37, not GRCh37.p13
    metadomain_dir = os.path.join(data_dir, "metadomains", args.genome_build.split(".")[0])

    for label, path in (("visualization", visualization_dir), ("metadomain", metadomain_dir)):
        if not os.path.isdir(path):
            raise SystemExit("{} directory not found: {}".format(label, path))

    paths = sorted(glob.glob(os.path.join(visualization_dir, "*", VISUALIZATION_FILE)))
    if args.limit:
        paths = paths[:args.limit]

    print("data dir     : {}".format(data_dir))
    print("genome build : {}".format(args.genome_build))
    print("transcripts  : {}".format(len(paths)))
    print("json parser  : {}".format(_JSON))
    if args.report:
        return
    if not paths:
        raise SystemExit("no transcripts found under {}".format(visualization_dir))

    started = time.time()
    print("loading ClinVar annotations from {} ...".format(metadomain_dir), flush=True)
    clinvar, clinvar_rows = load_clinvar(metadomain_dir)
    print("  {} records across {} (domain, consensus) keys in {:.1f}s".format(
        clinvar_rows, len(clinvar), time.time() - started), flush=True)
    if not clinvar:
        print("WARNING: no ClinVar records loaded - every accession list will be empty.",
              file=sys.stderr)

    shard_dir = args.out + ".shards"
    os.makedirs(shard_dir, exist_ok=True)
    jobs = [
        (index, bucket, os.path.join(shard_dir, "shard-{:04d}.tsv.gz".format(index)))
        for index, bucket in enumerate(chunk(paths, max(1, args.workers)))
    ]

    total = mismatches = 0
    try:
        with Pool(processes=len(jobs), initializer=_init_worker,
                  initargs=(clinvar, args.check)) as pool:
            for index, written, bad in pool.imap_unordered(_write_shard, jobs):
                total += written
                mismatches += bad
                print("  shard {:04d}: {} rows".format(index, written), flush=True)
    except ConsistencyError as exc:
        shutil.rmtree(shard_dir, ignore_errors=True)
        raise SystemExit("consistency check failed: {}".format(exc))

    # gzip members concatenate into a valid gzip file, so shards merge at disk
    # speed with no decompress/recompress step.
    print("merging {} shards -> {}".format(len(jobs), args.out), flush=True)
    with gzip.open(args.out, "wt", newline="", encoding="utf-8", compresslevel=6) as out:
        csv.writer(out, delimiter="\t", lineterminator="\n").writerow(COLUMNS)
    with open(args.out, "ab") as out:
        for _, _, shard_path in jobs:
            with open(shard_path, "rb") as shard:
                shutil.copyfileobj(shard, out, length=1024 * 1024)
    shutil.rmtree(shard_dir, ignore_errors=True)

    elapsed = time.time() - started
    print("\nwrote {} rows to {} in {:.1f}s ({:.0f} rows/s)".format(
        total, args.out, elapsed, total / elapsed if elapsed else 0))
    if mismatches:
        print("WARNING: {} placements where the accession lists disagreed with the "
              "prebuilt counts. Re-run with --check to stop at the first one.".format(mismatches))
    else:
        print("consistency: every placement's accession lists matched the prebuilt counts.")


if __name__ == "__main__":
    main()