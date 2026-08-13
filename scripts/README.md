# Scripts

Side-scripts for MetaDome operations, data processing and releases.

Scripts are grouped by what they need in order to run, since that determines
where they can be run.

| Directory                                    | Needs                          | Runs in                          |
| -------------------------------------------- | ------------------------------ | -------------------------------- |
| [`variant_annotation/`](variant_annotation/) | Flask app context + PostgreSQL | the `app` service |
| [`data_release/`](data_release/)             | prebuild output on disk only   | host or container, no database   |
| [`variant_analysis/`](variant_analysis/)     | a published release + bcftools/bedtools | host, no database      |

## Why the split

`variant_annotation/` runs inside the `app` service because it uses the Flask
application context and the repository layer, so it needs the database.
  
`data_release/` and `variant_analysis/` import nothing from `metadome` and read
only files from disk, so they run natively on the host. `data_release/` needs
Python; `variant_analysis/` needs `bcftools` and `bedtools`.
  
## Paths

Scripts resolve the data directory from `METADOME_DIR` in the repo-root `.env`,
the same variable `docker-compose.yml` substitutes into the `data` volume, so
one invocation works on the host and inside a container alike. `--data-dir`
overrides it; the fallback is `/usr/data`, the in-container mount point.

`.env` is committed with the key names. The values are machine-specific, so
keep your edits out of commits.

## Adding a script

Place it in the directory matching what it needs and add a `README.md` section
describing its inputs, outputs and invocation. If it needs an environment the
others do not, add it to the table above.
