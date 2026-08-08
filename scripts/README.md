# Scripts

Side-scripts for MetaDome operations, data processing and releases.

Scripts are grouped by what they need in order to run, since that determines
where they can be run.

| Directory                                    | Needs                          | Runs in                            |
| -------------------------------------------- | ------------------------------ | ---------------------------------- |
| [`variant_annotation/`](variant_annotation/) | Flask app context + PostgreSQL | the `app` service (`linux/amd64`)  |
| [`data_release/`](data_release/)             | prebuild output on disk only   | host or container, no database     |

## Why the split

`Dockerfile` pins the app image to `--platform=linux/amd64`, because it installs
x86-only BLAST+ and HMMER binaries. On Apple Silicon that image is emulated,
which costs time on anything CPU-bound. The database, redis and rabbitmq
services are `linux/arm64` and run natively.

Scripts under `data_release/` use neither BLAST nor HMMER nor the database, and
import nothing from `metadome`, so they need no Flask, SQLAlchemy or pysam and
run natively on the host.

## Paths

Scripts resolve the data directory from `METADOME_DIR` in the repo-root `.env`,
the same variable `docker-compose.yml` substitutes into the `data` volume, so
one invocation works on the host and inside a container alike. `--data-dir`
overrides it; the fallback is `/usr/data`, the in-container mount point.

`.env` is not committed. Copy `.env.example` and fill it in.

## Adding a script

Place it in the directory matching what it needs and add a `README.md` section
describing its inputs, outputs and invocation. If it needs an environment the
others do not, add it to the table above.
