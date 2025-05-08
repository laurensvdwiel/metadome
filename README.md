# Metadomain building pipeline

## Pre-requisites

### Required Files

This pipeline was initially built using:
**From [Gencode v45](https://www.gencodegenes.org/human/release_45.html):**
- `gencode.v45.pc_transcripts.fa`
- `gencode.v45.pc_translations.fa`
- `GRCh38.p14.genome.fa`
- `gencode.v45.annotation.gtf.gz`

**From [SwissProt 2025_01](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2025_01/):**
- `uniprot_sprot.dat`
- `uniprot_sprot_varsplic.fasta.gz`
- `uniprot_sprot.fasta.gz`

**From [PFAM v37.2](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.2/):**
- `Pfam-A.hmm`
- `Pfam-A.hmm.dat`


But you can update the download links in `metadomain_configs.yaml` and the name of the output file by modifying:  

```yaml
download_links:
  gencode:
    transcripts: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.pc_transcripts.fa.gz"
    translations: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.pc_translations.fa.gz"
    gtf: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz"
    genome: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.p14.genome.fa.gz"
    gencode_refseq: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.metadata.RefSeq.gz"
  uniprot:
    isoforms_fasta: "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz"
    swissprot_fasta: "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
  pfam:
    hmm: "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.4/Pfam-A.hmm.gz"
    hmm_dat: "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.4/Pfam-A.hmm.dat.gz"

step7_finaloutput:
    final_output: "data/results/step7_final_output/Gencodev47-PFAM37.4-Uniprot2025_metapositions.tsv"

```


---

### Required Software

_This pipeline relies on a number of tools beyond Python:_
- [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
- [HMMER](http://hmmer.org/)
- [Seqkit](https://bioinf.shenwei.me/seqkit/)
- [pfam_scan.pl](https://pfam-docs.readthedocs.io/en/latest/ftp-site.html)
- Python 

For easy sharing, they are currently included within a [pixi](https://prefix.dev/) environment. To obtain in and configure it correctly, clone this repo, then run:

```bash
git clone git@github.com:f-ferraro/pLMetadomains.git
cd pLMetadomains
pixi install 
```
---
---

## Running the pipeline 

You can then run the pipeline with the example command in `sbatch-metadomains.sh`
```bash
This scripts builds and annotate protein MetaDomains starting from Gencode, Uniprot, and PFAM input data.

Usage : sbatch-metadomains.sh -s <pixi.toml> -c <config> -r <clean-up> -t <int> -g <memory>
           -s Path to directory containing the manifest, the script, and where the analysis will be run.
           -c config.yaml that specifies inputs URLS to download the data from
           -r Boolean to specify to remove intermediate files at the end (default=False, alternative True)
```

Otherwise you can run directly the python script 
```bash
pixi run python metadomain_pipeline.py  -h 
usage: metadomain_pipeline.py [-h] --config CONFIG --cores CORES --working_dir_path WORKING_DIR_PATH

Complete pipeline to identify protein domains and create annotation files

options:
  -h, --help            show this help message and exit
  --config CONFIG       Path to yaml file containing information about files to download and save.
  --cores CORES         Number of available cores
  --working_dir_path WORKING_DIR_PATH
                        Path to directory containing the pixi.toml, this script and where analysis will be stored
```

