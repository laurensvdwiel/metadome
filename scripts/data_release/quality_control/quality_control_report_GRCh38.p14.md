# MetaDome data release quality control

Verification of the MetaDome GRCh38.p14 data release. Each check below reads only the files listed in the table and needs no database, prebuild output or running application, so it can be reproduced from the Zenodo record alone.

Release record: <https://zenodo.org/records/19376150>
Previous release: <https://zenodo.org/records/6625251>

| input | file |
| --- | --- |
| final dataset | `MetaDome_v2.0_GRCh38.p14_GENCODE-v45_UniProt-2025-01_Pfam-37.4_gnomAD-v4.1_ClinVar-2025-10-06_final-dataset-sw10.tsv.gz` |
| tolerance track | `MetaDome_v2.0_GRCh38.p14_GENCODE-v45_UniProt-2025-01_Pfam-37.4_gnomAD-v4.1_ClinVar-2025-10-06_derived-track-tolerance-landscape-sw10.bed.gz` |
| cross-build dataset | `MetaDome_v2.0_GRCh37.p13_GENCODE-v19_UniProt-2025-01_Pfam-37.4_gnomAD-r2.0.2_ClinVar-2025-10-06_final-dataset-sw10.tsv.gz` |

## GRCh38.p14 aggregation statistics

| statistic | value | share of codons |
| --- | --- | --- |
| codon and transcript combinations | 25213357 | |
| unique genomic codons | 11301412 | |
| covered by more than one gene | 68016 | 0.60% |
| transcripts disagree on window size | 0 | 0.00% |
| transcripts disagree on coverage | 37312 | 0.33% |
| transcripts disagree on dn/ds | 208327 | 1.84% |

The tolerance track reproduces from the final dataset: **yes** (11301412 codons compared, 11301412 identical).

## GRCh38.p14 against GRCh37.p13

Residues are matched by UniProt accession and position, which is independent of the genome assembly. The two builds annotate different GENCODE and gnomAD releases, so the final row is expected to be small: identical tolerance requires identical variant content in the sliding window, which different gnomAD releases rarely produce. It is not a measure of disagreement between the builds.

| statistic | value |
| --- | --- |
| residues in GRCh38.p14 | 19199448 |
| residues in GRCh37.p13 | 19212759 |
| shared | 17837395 (86.70% of all residues in either build) |
| GRCh38.p14 only | 1362053 |
| GRCh37.p13 only | 1375364 |
| shared residues with dn/ds equal to within 1e-09 | 386399 (2.17%) |

