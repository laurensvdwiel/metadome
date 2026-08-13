# MetaDome data release quality control

Verification of the MetaDome GRCh37.p13 data release. Each check below reads only the files listed in the table and needs no database, prebuild output or running application, so it can be reproduced from the Zenodo record alone.

Release record: <https://zenodo.org/records/19376150>
Previous release: <https://zenodo.org/records/6625251>

| input | file |
| --- | --- |
| final dataset | `MetaDome_v2.0_GRCh37.p13_GENCODE-v19_UniProt-2025-01_Pfam-37.4_gnomAD-r2.0.2_ClinVar-2025-10-06_final-dataset-sw10.tsv.gz` |
| tolerance track | `MetaDome_v2.0_GRCh37.p13_GENCODE-v19_UniProt-2025-01_Pfam-37.4_gnomAD-r2.0.2_ClinVar-2025-10-06_derived-track-tolerance-landscape-sw10.bed.gz` |
| v1.0.1 table | `metadome_data_full_n_transcripts_41772_20220508-011127.tsv.gz` |
| cross-build dataset | `MetaDome_v2.0_GRCh38.p14_GENCODE-v45_UniProt-2025-01_Pfam-37.4_gnomAD-v4.1_ClinVar-2025-10-06_final-dataset-sw10.tsv.gz` |

## GRCh37.p13 aggregation statistics

| statistic | value | share of codons |
| --- | --- | --- |
| codon and transcript combinations | 23768914 | |
| unique genomic codons | 10848533 | |
| covered by more than one gene | 34867 | 0.32% |
| transcripts disagree on window size | 0 | 0.00% |
| transcripts disagree on coverage | 41387 | 0.38% |
| transcripts disagree on dn/ds | 218254 | 2.01% |

The tolerance track reproduces from the final dataset: **yes** (10848533 codons compared, 10848533 identical).

## Comparison with MetaDome v1.0.1 (2022)

| statistic | value |
| --- | --- |
| codons in v1.0.1 | 10882508 |
| codons in this release | 10848533 |
| shared | 10605455 (95.32% of all codons in either release) |
| v1.0.1 only | 277053 |
| this release only | 243078 |
| shared codons with dn/ds equal to within 1e-09 | 10592889 (99.88%) |

## GRCh37.p13 against GRCh38.p14

Residues are matched by UniProt accession and position, which is independent of the genome assembly. The two builds annotate different GENCODE and gnomAD releases, so the final row is expected to be small: identical tolerance requires identical variant content in the sliding window, which different gnomAD releases rarely produce. It is not a measure of disagreement between the builds.

| statistic | value |
| --- | --- |
| residues in GRCh37.p13 | 19212759 |
| residues in GRCh38.p14 | 19199448 |
| shared | 17837395 (86.70% of all residues in either build) |
| GRCh37.p13 only | 1375364 |
| GRCh38.p14 only | 1362053 |
| shared residues with dn/ds equal to within 1e-09 | 386399 (2.17%) |

