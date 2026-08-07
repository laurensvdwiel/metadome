# Data release scripts

Build the public MetaDome data release: the final dataset and its derived
tracks, one set per genome build.

These scripts read only prebuild output from disk. They use no database, no
running application, and import nothing from `metadome`, so they run on the
host without the `app` image.

`orjson` is optional and makes JSON parsing several times faster. Without it
the scripts use the standard library.

## Pipeline

```
prebuild output on disk
        │
        ▼  generate_final_dataset.py        (per genome build)
   ..._final-dataset-sw10.tsv.gz
        │
        ▼  generate_derived_tracks.py
   ..._derived-track-tolerance-landscape-sw10.bed
   ..._derived-track-metadomain-clinvar.bed
   ..._derived-track-pfam-domain-coverage.bed
```

Every track is built from the final dataset, so they always describe the same
data.

## `generate_final_dataset.py`

Joins two on-disk sources per genome build:

- `metadome_visualization/<build>/<transcript>/metadome_visualization.json` —
  codon position, strand, protein coordinates, tolerance, per-domain counts
- `metadomains/<build>/<PFxxxxx>/metadomain_snv_annotation_clinvar` — ClinVar
  accessions, keyed by domain consensus position

Output is 30 columns, one row per transcript, codon fragment and domain
placement. A codon outside every Pfam domain gets a row with the domain columns
empty, and a codon split across an exon boundary yields one row per fragment.

Genomic positions are 1-based inclusive. Protein positions are 1-based from the
initiator methionine. The two sources are joined on
`(domain_id, consensus_pos)`.

Homologous counts and accessions describe the evolutionarily equivalent
positions in other genes; ClinVar records on the annotated codon itself are
reported separately in `clinvar_P_at_position` and `clinvar_LP_at_position`.
Counts are of distinct single nucleotide variants, so an accession list can be
shorter than its count when one ClinVar record appears in more than one reading
frame.

```bash
# pilot first, to validate the join on a small slice
python scripts/data_release/generate_final_dataset.py \
    --genome-build GRCh37.p13 \
    --out /path/to/pilot.tsv.gz \
    --limit 500 --workers 4 --check

# full run
python scripts/data_release/generate_final_dataset.py \
    --genome-build GRCh37.p13 \
    --out /path/to/MetaDome_v2.0_GRCh37.p13_..._final-dataset-sw10.tsv.gz \
    --workers 4 --check
```

`--data-dir` overrides `METADOME_DIR` from the repo-root `.env`; the fallback is
`/usr/data`. Keep `--workers` below the core count while a prebuild is running.

### `--check`

The prebuild writes `pathogenic_variant_count_per_clinsig` into each
visualization JSON. This script derives the accession lists independently from
the metadomain CSVs. `--check` verifies that the counts summed across a
residue's consensus positions match the prebuilt totals, and stops at the first
disagreement.

Two independent paths to the same number, so a fault in the consensus-position
join or in the self-exclusion rule surfaces as a failure rather than as a
plausible-looking file. The full run reports a count either way.

## `generate_derived_tracks.py`

```bash
python scripts/data_release/generate_derived_tracks.py \
    --final-dataset /path/to/..._final-dataset-sw10.tsv.gz \
    --genome-build GRCh37.p13 \
    --out-dir /path/to/release \
    --prefix MetaDome_v2.0_GRCh37.p13_GENCODE-v19_UniProt-2025-01_Pfam-37.4_gnomAD-r2.0.2_ClinVar-2025-10-06
```

Emits three BED files and prints the aggregation statistics. Each file opens
with a commented header carrying the release stem, the assembly and the Zenodo
record.

The tolerance track holds one row per genomic codon. Where transcripts overlap
a codon and their scores disagree, the median is used; the printed statistics
report how often that happens.

`--skip-tolerance` and `--skip-domains` build a subset.

## `metadomain_clinvar.as` and `pfam_coverage.as`

autoSql schemas describing the extra columns of the two domain tracks,
`bed9+11` and `bed6+5`. BED defines names only for its first twelve columns, so
these files document the rest: each extra field's name, type and meaning.

They also validate the tracks. The commented header is stripped first:

```bash
grep -v '^#' <prefix>_derived-track-metadomain-clinvar.bed \
  | bedToBigBed -type=bed9+11 -as=metadomain_clinvar.as -tab \
      stdin <chrom.sizes> <prefix>_derived-track-metadomain-clinvar.bb

grep -v '^#' <prefix>_derived-track-pfam-domain-coverage.bed \
  | bedToBigBed -type=bed6+5 -as=pfam_coverage.as -tab \
      stdin <chrom.sizes> <prefix>_derived-track-pfam-domain-coverage.bb
```

## Counting

Two keys identify a position, and they answer different questions.

`(chrom, pos_start, pos_stop, strand)` identifies a genomic interval. Use it for
totals: overlapping genes place two protein positions over one set of variants,
and only the genomic key collapses them.

`(uniprot_ac, protein_pos, domain_id, consensus_pos)` identifies a residue and
is independent of the assembly, so it compares positions between genome builds.
It also collapses a split codon, which the genomic key counts as two fragments.

Count columns are per row. Summing them across the dataset multiplies each
position by the number of transcripts covering it.
