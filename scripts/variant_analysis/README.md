# Variant analysis scripts
  
Analyses that annotate external variant sets against a published MetaDome data
release.
  
These read only files from the release, so they need no database, no prebuild
output and no running application. Anyone holding the Zenodo record can
reproduce the result with `bcftools` and `bedtools`.
  
## `clinvar_analysis.sh`
  
Finds ClinVar variants of uncertain significance that sit at positions where
homologous Pfam domains carry pathogenic missense variation, i.e. positions
where evidence can be transferred from another genomic position through homologous protein domains.
  
### Inputs
  
| Input | Source |
| --- | --- |
| `--clinvar` | ClinVar VCF for the assembly being analysed |
| `--release-dir` + `--prefix` | `_derived-track-pfam-domain-coverage.bed.gz`, `_derived-track-metadomain-clinvar.bed.gz` and `_final-dataset-sw10.tsv.gz` from the release |
  
The ClinVar VCF must match the assembly of the release; the two are joined on
genomic coordinates and nothing detects a mismatch.
  
### Running
  
```bash
./clinvar_analysis.sh \
    --release-dir  /path/to/release \
    --prefix       MetaDome_v2.0_GRCh38.p14_GENCODE-v45_UniProt-2025-01_Pfam-37.4_gnomAD-v4.1_ClinVar-2025-10-06 \
    --clinvar      /path/to/ClinVar/GRCh38/clinvar_20251006.vcf.gz \
    --genome-build GRCh38.p14 \
    --out-dir      /path/to/output
```

`--clnsig` selects a different ClinVar significance class, so the same script
produces sensitivity arms over `Pathogenic/Likely_pathogenic` or
`Conflicting_classifications_of_pathogenicity`. `--keep-intermediates` retains
the BED files and raw intersects.
  
### Output
  
`vus_metadomain_<build>.tsv` holds one row per variant per domain placement,
restricted to variants with at least one homologous pathogenic or likely
pathogenic missense variant. Columns are the variant as ClinVar reports it,
then the MetaDome annotation: gene symbol, UniProt accession and position, Pfam
accession and consensus position, the homologous counts, the accessions behind
those counts, and a link to the position in MetaDome.
  
A variant covered by several domain placements or isoforms appears more than
once, so per-variant figures come from distinct `clinvar_id`, not row counts.
  
`vus_metadomain_<build>_summary.tsv` records the counts the table cannot: how
many variants were selected, and how many fall inside a Pfam domain. Variants
without homologous evidence are excluded from the table, so its denominators
live here.
  
### Coordinates
  
The release tracks are 1-based inclusive on both ends: `chr1 878710 878712` is
the three bases `878710` to `878712`. `bedtools` reads column 2 as 0-based
half-open, so the script decrements it before intersecting. Without that shift
every variant on a codon's first base is silently missed.
  
The script verifies this held. A codon's own ClinVar records are excluded from
its homologous counts when the release is built, so a variant must never appear
among its own accessions. The reported `self-contamination` count must be zero;
anything else means the coordinate conversion is wrong and the output should be
discarded.

## `plot_clinvar_analysis.py`
  
Draws a plotted figure from the output of `clinvar_analysis.sh` for one genome build. 

### Running
  
`pandas`, `seaborn` and `matplotlib` are not application dependencies and are
not in `requirements.txt`. `Dockerfile_plots` provides them:
  
```bash
docker build -f Dockerfile_plots -t metadome-plots .
  
docker run --rm -u "$(id -u):$(id -g)" \
    -v /path/to/output:/work -v "$PWD:/scripts:ro" metadome-plots \
    python /scripts/plot_clinvar_analysis.py \
        --table vus_metadomain_GRCh38.p14.tsv \
        --summary vus_metadomain_GRCh38.p14_summary.tsv \
        --genome-build GRCh38.p14 \
        --counts
```

`-u` keeps the output owned by the invoking user rather than root. `/work` is
the directory holding the tables, so `--table` and `--summary` are relative to
it and the figures are written beside them.
  
### Output
  
`vus_metadomain_<build>.pdf` and `.png` at 300 dpi. `--counts` additionally
writes `vus_metadomain_<build>_panelB_counts.tsv`, the raw counts behind panel
B, which partition the candidate set.
  
The thresholds printed on completion are the check: they are computed by pandas
and must match the counts `clinvar_analysis.sh` reports.
