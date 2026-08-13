# MetaDome
  
MetaDome analyses genetic variants by aggregating homologous human protein
domains. Variation is projected onto Pfam domain consensus positions, so
variants in different genes that sit on evolutionarily equivalent residues can
be compared with one another.
  
The public platform runs at [www.metadome.app](https://www.metadome.app/metadome).
This repository is everything needed to run your own instance.

## Repository layout
  
| Path | Contains |
| --- | --- |
| `metadome/` | The application: Flask web layer, Celery tasks, domain model, repositories |
| `metadome_prebuild_GRCh37/`, `metadome_prebuild_GRCh38/` | Standalone pixi pipelines that build a mapping database for one genome build |
| `scripts/data_release/` | Builds the public data release and its genome-browser tracks |
| `scripts/data_release/quality_control/` | Verifies a release and writes the reports committed beside it |
| `scripts/variant_annotation/` | Annotates externally supplied variants with meta-domain positions |
| `tests/` | Unit and integration tests |

## Requirements

Docker [here](https://www.docker.com/get-docker) and Docker-compose [here](https://docs.docker.com/compose/install/#install-compose). Everything else runs in containers, and the images build natively on both x86-64 and arm64.

## Configuration

`.env` holds the four directories the containers mount. Set each to an absolute 
path on your machine:
  
    METADOME_DIR=/absolute/path/to/metadome
    CLINVAR_DIR=/absolute/path/to/ClinVar
    GNOMAD_DIR=/absolute/path/to/gnomAD
    METADOME_POSTGRES_DB_DIR=/absolute/path/to/an/empty/directory

`METADOME_POSTGRES_DB_DIR` must be empty the first time the database starts.
`METADOME_DIR` is where prebuilt output is read from and written to.
  
`.env` is tracked so the key names travel with the repository. The values are
machine-specific, so keep your edits out of commits.

### Credentials
  
Before exposing an instance publicly, replace the defaults:
  
- `POSTGRES_PASSWORD`, `POSTGRES_USER`, `POSTGRES_DB` in `metadome/database.env`
- `SECRET_KEY_CRED` in `metadome/flask_app_credentials.py`
- the relay host settings in `metadome/mailserver.env`

## Reference data 
  
The application reads two variant sources at runtime:
  
| Path | File                                                   |
| --- |--------------------------------------------------------|
| `$CLINVAR_DIR/GRCh37/` | `clinvar_YYYYMMDD.vcf.gz` and its `.tbi`               |
| `$CLINVAR_DIR/GRCh38/` | `clinvar_YYYYMMDD.vcf.gz` and its `.tbi`               |
| `$GNOMAD_DIR/GRCh37/` | `gnomad.exomes.V2XXX.sites.vcf.gz` and its `.tbi`      |
| `$GNOMAD_DIR/GRCh38/` | `gnomad.joint.v4.X.sites.exomes.vcf.gz` and its `.tbi` |
  
Filenames are set in `metadome/default_settings.py`; change them there to use
other releases. Only variants passing gnomAD's `PASS` filter are counted.

#### ClinVar Version support

MetaDome was tested with version 20251006 for GRCh37 and GRCh38, both [available here](ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/), but later versions may work as well.


#### gnomAD Version support

MetaDome was tested with version v2 for GRCh37 and v4.1 (joint exomes) for GRCh38, both  [available here](https://gnomad.broadinstitute.org/data), but later versions may work as well as they come available.

#### GENCODE, UniProt and Pfam 

GENCODE, UniProt and Pfam are not read by the application. They are inputs to
the prebuild pipelines, which fetch their own copies.


## Setting up
  
### 1. Obtain a mapping database
  
The mapping database links every coding genomic position to a transcript, a
UniProt residue and, where the residue falls inside one, a Pfam domain
consensus position. Everything else is built on it.
  
**Download the current release.** Both genome builds are published at [Zenodo record 19376150](https://zenodo.org/records/19376150):
  
    MetaDome_v2.0_GRCh37.p13_GENCODE-v19_UniProt-2025-01_Pfam-37.4_prebuild-mapping-database.csv.gz
    MetaDome_v2.0_GRCh38.p14_GENCODE-v45_UniProt-2025-01_Pfam-37.4_prebuild-mapping-database.csv.gz
  
**Or build your own.** See  
[`metadome_prebuild_GRCh37/README.md`](metadome_prebuild_GRCh37/README.md) and 
[`metadome_prebuild_GRCh38/README.md`](metadome_prebuild_GRCh38/README.md).

### 2. Load it into the database
  
`populate_database.sh` brings up a temporary stack, loads the CSV through
`metadome.batch_load`, and tears the stack down again:
  
```bash
./populate_database.sh \
    --data-dir /path/to/the/directory/holding/the/csv \
    --file MetaDome_v2.0_GRCh38.p14_GENCODE-v45_UniProt-2025-01_Pfam-37.4_prebuild-mapping-database.csv.gz
```
  
`--file` is relative to `--data-dir`, which is mounted into the container at
`/prebuild_data`. Run it once per genome build; both coexist in one database.
  
The mapping table holds tens of millions of rows per build, so this is a long
load and needs disk headroom for the database plus its indexes.

### 3. Prebuild visualisations and meta-domains (optional)
  
MetaDome computes a transcript's tolerance landscape and meta-domain annotation
on demand, but the first request for a transcript then takes minutes.
Prebuilding stores the result so requests are served immediately. A deployment
works without this step; it is only slower.
  
**Download the current release.** From the same Zenodo record, per genome build:

    MetaDome_v2.0_<build>_..._prebuild-visualizations.tar.zst
    MetaDome_v2.0_<build>_..._prebuild-metadomains.tar.zst

The two archives have different roots, so they extract to different places
under `METADOME_DIR`:

```bash
# roots at GRCh37.p13/ and GRCh38.p14/
zstd -dc MetaDome_v2.0_..._prebuild-visualizations.tar.zst \
  | tar -xf - -C "$METADOME_DIR/metadome_visualization/"
  
# roots at GRCh37/ and GRCh38/
zstd -dc MetaDome_v2.0_..._prebuild-metadomains.tar.zst \
  | tar -xf - -C "$METADOME_DIR/metadomains/"
```

The release these archives are cut from is verified separately, and the reports are committed alongside the code that produces them:
  
- [`quality_control_report_GRCh37.p13.md`](scripts/data_release/quality_control/quality_control_report_GRCh37.p13.md)
- [`quality_control_report_GRCh38.p14.md`](scripts/data_release/quality_control/quality_control_report_GRCh38.p14.md)
  
Each records the aggregation statistics recomputed from the released dataset,
confirms the tolerance track reproduces from it codon for codon, and compares
the release against MetaDome v1.0.1 (2022) and against the other genome build.

**Or build them yourself**, once the database is populated:

```bash
./prebuild_all_after_populate_database.sh --workers 8
```
  
Options are forwarded to `metadome.prebuild_all`: 
  
| Option | Effect |
| --- | --- | 
| `--report` | Audit what is pending and build nothing |
| `--workers N` | Parallel worker processes |
| `--limit N` | Build at most N pending items, then stop |
| `--genome-build B` | Restrict to one build; repeatable |
| `--overwrite` | Rebuild everything, including what already exists |
| `--stale-days D` | Rebuild entries older than D days |
  
This reuses a running stack and needs only the database service, so it is safe
against a live instance. A full build across both genome builds takes days.
  
## Running the platform 
  
```bash
docker compose up -d
```
  
The application is served by gunicorn on port 5000. Stop it with `docker compose stop`.
  
## Tests
  
```bash
./run_tests.sh                              # unit tests 
./run_tests.sh --tests tests/integration    # integration tests
```

## Data release and quality control
  
`scripts/data_release/` builds the published dataset and three genome-browser
tracks per genome build from prebuild output alone, with no database and no
running application. `scripts/data_release/quality_control/` verifies a release
against itself and against the 2022 release. Both have their own READMEs.

## Citing MetaDome

Please cite the MetaDome manuscript :

```
MetaDome: Pathogenicity analysis of genetic variants
 through aggregation of homologous human protein domains.
Laurens Wiel, Coos Baakman, Daan Gilissen, Joris A. Veltman,
 Gert Vriend and Christian Gilissen.
Human Mutation. (2019) 1-9, 10.1002/humu.23798
```

## Contact

If you want to provide feedback please have a look at our
[existing issues][1] (and if necessary, create a new issue).

[1]: https://github.com/laurensvdwiel/metadome/issues
