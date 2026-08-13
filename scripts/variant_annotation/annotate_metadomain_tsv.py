import argparse
import csv
import sys

from metadome.factory import create_app
from metadome.domain.services.external.meta_domain_position_annotation import annotate_variants_with_metadomain

DEFAULT_BATCH_SIZE = 25


def parse_variant(variant_str):
    """Parse 'chr1:123456_A_G' into (chr, pos, ref)."""
    chrom, rest = variant_str.split(':', 1)
    pos_str, ref, _alt = rest.split('_', 2)
    return chrom, int(pos_str), ref


def main():
    parser = argparse.ArgumentParser(
        description="Annotate a TSV file with MetaDomain positions directly via app context."
    )
    parser.add_argument("input", help="Input TSV file")
    parser.add_argument("output", help="Output TSV file")
    parser.add_argument(
        "--genome-build",
        required=True,
        help="Genome build to use for filtering, e.g. GRCh38 or GRCh38.p14",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=DEFAULT_BATCH_SIZE,
        help=f"Number of variants per processing batch (default: {DEFAULT_BATCH_SIZE})",
    )
    parser.add_argument(
        "--max-rows",
        type=int,
        default=None,
        help="Only process and write the first N rows",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print progress information",
    )
    parser.add_argument(
        "--stop-on-error",
        action="store_true",
        help="Stop immediately on the first failed batch",
    )
    args = parser.parse_args()

    with open(args.input, "r", newline="") as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        if "variant" not in fieldnames or "gene_id" not in fieldnames:
            raise RuntimeError("Input TSV must contain 'variant' and 'gene_id' columns")
        rows = list(reader)

    if args.max_rows is not None:
        rows = rows[:args.max_rows]

    for extra_field in ["MetaDomainPositions", "MetaDomainStatus", "RefMatchStatus"]:
        if extra_field not in fieldnames:
            fieldnames.append(extra_field)

    parsed_rows = []
    for index, row in enumerate(rows, start=1):
        try:
            chrom, pos, ref = parse_variant(row["variant"])
            parsed_rows.append({
                "row_index": index - 1,
                "chr": chrom,
                "pos": pos,
                "ref": ref,
                "gene_id": row["gene_id"],
                "genome_build": args.genome_build,
            })
        except Exception as exc:
            if args.verbose:
                print(
                    f"[WARN] Could not parse row {index}: variant={row.get('variant')} "
                    f"gene_id={row.get('gene_id')} error={exc}",
                    file=sys.stderr,
                    flush=True,
                )
            parsed_rows.append(None)

    annotations = [
        {
            "MetaDomainPositions": "",
            "MetaDomainStatus": "",
            "RefMatchStatus": "",
        }
        for _ in rows
    ]
    valid_entries = [entry for entry in parsed_rows if entry is not None]

    if args.verbose:
        print(
            f"[INFO] Loaded {len(rows)} rows, {len(valid_entries)} parseable variants, "
            f"batch_size={args.batch_size}, genome_build={args.genome_build}",
            file=sys.stderr,
            flush=True,
        )

    app = create_app()

    with app.app_context():
        total_batches = (len(valid_entries) + args.batch_size - 1) // args.batch_size

        for batch_number, start in enumerate(range(0, len(valid_entries), args.batch_size), start=1):
            batch = valid_entries[start:start + args.batch_size]

            try:
                if args.verbose:
                    print(
                        f"[INFO] Processing batch {batch_number}/{total_batches} "
                        f"with {len(batch)} variants",
                        file=sys.stderr,
                        flush=True,
                    )

                results = annotate_variants_with_metadomain(batch)

                result_map = {
                    (
                        result["chr"],
                        int(result["pos"]),
                        result.get("ref"),
                        result["gene_id"].split(".", 1)[0] if result.get("gene_id") else result.get("gene_id"),
                        result["genome_build"],
                    ): {
                        "MetaDomainPositions": result.get("MetaDomainPositions", ""),
                        "MetaDomainStatus": result.get("MetaDomainStatus", ""),
                        "RefMatchStatus": result.get("RefMatchStatus", ""),
                    }
                    for result in results
                }

                for entry in batch:
                    key = (
                        entry["chr"],
                        int(entry["pos"]),
                        entry.get("ref"),
                        entry["gene_id"].split(".", 1)[0] if entry.get("gene_id") else entry.get("gene_id"),
                        entry["genome_build"],
                    )
                    annotations[entry["row_index"]] = result_map.get(
                        key,
                        {
                            "MetaDomainPositions": "",
                            "MetaDomainStatus": "no_mapping",
                            "RefMatchStatus": "not_checked",
                        },
                    )

            except Exception as exc:
                print(
                    f"[ERROR] Failed batch {batch_number}/{total_batches}: {exc}",
                    file=sys.stderr,
                    flush=True,
                )
                if args.stop_on_error:
                    raise

            if args.max_rows is not None and (start + len(batch)) >= args.max_rows:
                break

    with open(args.output, "w", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for row, annotation in zip(rows, annotations):
            row["MetaDomainPositions"] = annotation["MetaDomainPositions"]
            row["MetaDomainStatus"] = annotation["MetaDomainStatus"]
            row["RefMatchStatus"] = annotation["RefMatchStatus"]
            writer.writerow(row)

    print(
        f"[INFO] Done. Wrote {len(rows)} rows to {args.output}",
        file=sys.stderr,
        flush=True,
    )


if __name__ == "__main__":
    main()