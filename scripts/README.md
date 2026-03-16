# Scripts

This directory contains internal utility scripts for MetaDome operations and data processing.

## `annotate_metadomain_tsv.py`

Annotates a TSV file with MetaDomain information using the MetaDome application context directly.

### What it does

The script expects an input TSV that contains at least these columns:

- `variant`
- `gene_id`

The `variant` column must use this format:
`chr1:123456_A_G`

From this, the script extracts:

- chromosome
- position
- reference allele

It then annotates each row with these output columns:

- `MetaDomainPositions`
- `MetaDomainStatus`
- `RefMatchStatus`

### Output columns

#### `MetaDomainPositions`

A semicolon-separated list of metadomain positions:
`PF00001:42` or `PF00002:15;PF00002:18`


#### `MetaDomainStatus`

One of:

- `metadomain_found`
- `mapping_no_metadomain`
- `no_mapping`
- `invalid_input`

#### `RefMatchStatus`

One of:

- `direct_match`
- `reverse_complement_match`
- `mismatch`
- `not_checked`

### Example command

Run the script through Docker Compose so it reuses the configured app environment and service network:
``` 
docker compose run --rm
-v "<HOST_DATA_DIR>":/data
app
python -m scripts.annotate_metadomain_tsv
/data/input_variants.tsv
/data/output_annotated.tsv
--genome-build GRCh38.p14
--batch-size 500
--verbose
--stop-on-error
```

### Arguments

### Notes

- This script is intended for internal/batch use.
- It does not rely on the public MetaDome API.
- It runs inside the Flask application context and uses the repository/service layer directly.
- The external data directory must be mounted explicitly, for example with: `-v "<HOST_DATA_DIR>":/data`
