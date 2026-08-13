# Variant annotation scripts

Scripts that annotate externally supplied variants using MetaDome's repository
and service layer.

These require the Flask application context and a running PostgreSQL. Run them
through `docker compose` so they reuse the configured app environment and
service network.

## `annotate_metadomain_tsv.py`

Annotates a TSV of variants with the Pfam meta-domain positions they fall on.

Given a variant list, it reports which variants sit at evolutionarily
equivalent positions across the human proteome, so variation scattered across
different genes can be compared at the same domain consensus position.

### Input

A TSV containing at least:

- `variant`, formatted `chr1:123456_A_G`
- `gene_id`

Chromosome, position and reference allele are parsed from `variant`.

### Output columns

| Column                | Meaning                                                                        |
| --------------------- | ------------------------------------------------------------------------------ |
| `MetaDomainPositions` | Semicolon-separated meta-domain positions, e.g. `PF00002:15;PF00002:18`        |
| `MetaDomainStatus`    | `metadomain_found`, `mapping_no_metadomain`, `no_mapping` or `invalid_input`   |
| `RefMatchStatus`      | `direct_match`, `reverse_complement_match`, `mismatch` or `not_checked`        |

Minus-strand genes store transcript-oriented codons, so a reference allele that
appears not to match may be the reverse complement. `RefMatchStatus` reports
which of the two occurred.

### Running

```bash
docker compose run --rm \
  -v "<HOST_DATA_DIR>":/data \
  app \
  python -m scripts.variant_annotation.annotate_metadomain_tsv \
    /data/input_variants.tsv \
    /data/output_annotated.tsv \
    --genome-build GRCh38.p14 \
    --batch-size 500 \
    --verbose \
    --stop-on-error
```

The external data directory is mounted explicitly with `-v`.

### Notes

- For internal and batch use; it does not go through the public MetaDome API.
- It runs inside the Flask application context and uses the repository and
  service layer directly, so the database must be running.
- The meta-domain positions it returns can be intersected against the
  pathogenic records in the release final dataset.
