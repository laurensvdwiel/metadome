#!/bin/bash
#
# Annotate ClinVar variants of uncertain significance with the pathogenic
# missense variation observed at evolutionarily equivalent positions in
# homologous Pfam domains, and with the evidence that qualifies or contradicts
# it: benign variation at the same meta-domain position, and pathogenic
# variation at the variant's own residue.
#
# Output is one row per variant per domain placement: a variant sitting in two
# overlapping domains appears twice. Count distinct clinvar_id for per-variant
# figures.

set -euo pipefail
  
RELEASE_DIR=""
PREFIX=""
CLINVAR_VCF=""
GENOME_BUILD=""
OUT_DIR="."
CLNSIG="Uncertain_significance"
# Benign classes used to detect conflicting evidence. ClinVar's combined
# Benign/Likely_benign is included, unlike its pathogenic counterpart: omitting
# it would miss conflicts and overstate the eligible set.
BENIGN_EXPR='INFO/CLNSIG="Benign" || INFO/CLNSIG="Likely_benign" || INFO/CLNSIG="Benign/Likely_benign"'
# Pathogenic classes, matching the release convention: Pathogenic and
# Likely_pathogenic exactly, excluding ClinVar's combined class.
PATHOGENIC_EXPR='INFO/CLNSIG="Pathogenic" || INFO/CLNSIG="Likely_pathogenic"'
KEEP_INTERMEDIATES=0
  
function show_usage {
    cat <<'USAGE'
Usage: clinvar_analysis.sh [options]
  
Required:
  -r, --release-dir DIR    Directory holding the release files
  -p, --prefix STEM        Release filename stem, e.g.
                           MetaDome_v2.0_GRCh38.p14_GENCODE-v45_UniProt-2025-01_Pfam-37.4_gnomAD-v4.1_ClinVar-2025-10-06
  -c, --clinvar FILE       ClinVar VCF for the same assembly (.vcf.gz)
  -b, --genome-build NAME  Assembly label used in the output filename, e.g. GRCh38.p14
  
Optional:
  -o, --out-dir DIR        Where to write output (default: current directory)
  -s, --clnsig VALUE       ClinVar significance of the query variants
                           (default: Uncertain_significance). The benign and
                           pathogenic comparison sets are fixed.
  -k, --keep-intermediates Retain the BED and intersect files
  -h, --help               Show this message
  
The release files are located as:
  <release-dir>/<prefix>_derived-track-pfam-domain-coverage.bed.gz
  <release-dir>/<prefix>_derived-track-metadomain-clinvar.bed.gz
  <release-dir>/<prefix>_final-dataset-sw10.tsv.gz
USAGE
}
  
while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--release-dir)       RELEASE_DIR="$2"; shift 2 ;;
        -p|--prefix)            PREFIX="$2"; shift 2 ;;
        -c|--clinvar)           CLINVAR_VCF="$2"; shift 2 ;;
        -b|--genome-build)      GENOME_BUILD="$2"; shift 2 ;;
        -o|--out-dir)           OUT_DIR="$2"; shift 2 ;;
        -s|--clnsig)            CLNSIG="$2"; shift 2 ;;
        -k|--keep-intermediates) KEEP_INTERMEDIATES=1; shift ;;
        -h|--help)              show_usage; exit 0 ;;
        *)                      echo "Unknown option: $1" >&2; show_usage; exit 1 ;;
    esac
done
  
for required in RELEASE_DIR PREFIX CLINVAR_VCF GENOME_BUILD; do
    if [ -z "${!required}" ]; then
        echo "ERROR: --${required,,} is required." | tr '_' '-' >&2
        show_usage
        exit 1
    fi
done
  
for tool in bcftools bedtools awk sort; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: $tool is not on the PATH." >&2; exit 1; }
done
  
COVERAGE="$RELEASE_DIR/${PREFIX}_derived-track-pfam-domain-coverage.bed.gz"
EVIDENCE="$RELEASE_DIR/${PREFIX}_derived-track-metadomain-clinvar.bed.gz"
DATASET="$RELEASE_DIR/${PREFIX}_final-dataset-sw10.tsv.gz"
  
for f in "$CLINVAR_VCF" "$COVERAGE" "$EVIDENCE" "$DATASET"; do
    [ -f "$f" ] || { echo "ERROR: missing input: $f" >&2; exit 1; }
done
  
mkdir -p "$OUT_DIR"
WORK="$OUT_DIR/.clinvar_analysis_${GENOME_BUILD}"
mkdir -p "$WORK"
FINAL="$OUT_DIR/vus_metadomain_${GENOME_BUILD}.tsv"

# The release tracks name chromosomes chr1..chrY. ClinVar usually does not, so
# the prefix is added when absent rather than assumed either way.
if gzip -dc "$CLINVAR_VCF" | grep -v '^#' | head -1 | cut -f1 | grep -q '^chr'; then
    CHR_PREFIX=""
else
    CHR_PREFIX="chr"
fi

echo "genome build : $GENOME_BUILD"
echo "clinvar      : $(basename "$CLINVAR_VCF")"
echo "significance : $CLNSIG"
echo "chrom prefix : ${CHR_PREFIX:-<none added>}"

# ---------------------------------------------------------------------------
# 1. The variants of interest, as 0-based half-open BED.
# ---------------------------------------------------------------------------
echo "selecting missense SNVs ..."
bcftools view -v snps -i "INFO/CLNSIG=\"${CLNSIG}\" && INFO/MC~\"missense_variant\"" "$CLINVAR_VCF" \
  | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNHGVS\n' \
  | awk -F'\t' -v p="$CHR_PREFIX" 'BEGIN{OFS="\t"}{print p$1, $2-1, $2, $3, $4, $5, $6, $7, $8}' \
  | sort -k1,1 -k2,2n > "$WORK/vus.bed"

# ---------------------------------------------------------------------------
# 2. Release tracks, decompressed and sorted.
#
# The tracks are 0-based half-open BED.
# ---------------------------------------------------------------------------
echo "preparing release tracks ..."
gzip -dc "$COVERAGE" | grep -v '^#' | sort -k1,1 -k2,2n > "$WORK/coverage.bed"
gzip -dc "$EVIDENCE" | grep -v '^#' | sort -k1,1 -k2,2n > "$WORK/evidence.bed"

# ---------------------------------------------------------------------------
# 3a. Two independent intersects.
#
# Coverage gives every variant inside a Pfam domain, which is the denominator.
# Evidence gives those with pathogenic missense variation at the equivalent
# position elsewhere. They are joined on the placement in step 5 rather than
# chained, because bedtools matches on coordinates alone and would otherwise
# pair one domain placement with another placement's evidence.
# ---------------------------------------------------------------------------
echo "intersecting with Pfam domain coverage ..."
bedtools intersect -a "$WORK/vus.bed" -b "$WORK/coverage.bed" -wa -wb -sorted > "$WORK/in_domain.tsv"
echo "intersecting with meta-domain ClinVar evidence ..."
bedtools intersect -a "$WORK/vus.bed" -b "$WORK/evidence.bed" -wa -wb -sorted > "$WORK/with_evidence.tsv"

# ---------------------------------------------------------------------------
# 3b. Benign variation at the same meta-domain positions.
#
# The calibration for PM5 applies to meta-positions whose pathogenic evidence
# is not contradicted: a position carrying both pathogenic and benign
# observations is conflicting. MetaDome annotates only pathogenic and likely
# pathogenic variation, so the benign side is taken from ClinVar directly and
# mapped onto meta-positions through the same Pfam coverage track.
# ---------------------------------------------------------------------------
echo "selecting benign missense SNVs ..."
bcftools view -v snps -i "(${BENIGN_EXPR}) && INFO/MC~\"missense_variant\"" "$CLINVAR_VCF" \
  | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNHGVS\n' \
  | awk -F'\t' -v p="$CHR_PREFIX" 'BEGIN{OFS="\t"}{print p$1, $2-1, $2, $3, $4, $5, $6, $7, $8}' \
  | sort -k1,1 -k2,2n > "$WORK/benign.bed"

echo "intersecting benign variants with Pfam domain coverage ..."
bedtools intersect -a "$WORK/benign.bed" -b "$WORK/coverage.bed" -wa -wb -sorted \
  > "$WORK/benign_in_domain.tsv"

# ---------------------------------------------------------------------------
# 3c. Pathogenic variation at the variant's own residue.
#
# A pathogenic missense variant at the same amino acid position in the same
# protein is stronger evidence than the paralogue route and is assessed first.
# Variants carrying it would qualify under PM5 regardless of meta-domain
# evidence, so isolating those without it identifies where the meta-domain adds
# evidence that is otherwise unavailable.
# ---------------------------------------------------------------------------
echo "selecting pathogenic missense SNVs ..."
bcftools view -v snps -i "(${PATHOGENIC_EXPR}) && INFO/MC~\"missense_variant\"" "$CLINVAR_VCF" \
  | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNHGVS\n' \
  | awk -F'\t' -v p="$CHR_PREFIX" 'BEGIN{OFS="\t"}{print p$1, $2-1, $2, $3, $4, $5, $6, $7, $8}' \
  | sort -k1,1 -k2,2n > "$WORK/pathogenic.bed"

echo "intersecting pathogenic variants with Pfam domain coverage ..."
bedtools intersect -a "$WORK/pathogenic.bed" -b "$WORK/coverage.bed" -wa -wb -sorted \
  > "$WORK/same_protein_in_domain.tsv"

# ---------------------------------------------------------------------------
# 4. UniProt accession to gene symbol.
#
# Neither track carries a gene symbol, so it comes from the final dataset. The
# symbol is taken from MetaDome rather than from ClinVar's GENEINFO, which is
# ClinVar's own gene assignment and can disagree with the mapping used here.
# GENCODE falls back to an ENSG identifier where no HGNC symbol exists, so a
# real symbol is preferred; accessions shared by several genes keep all of
# them, joined with "|", rather than one being chosen silently.
# ---------------------------------------------------------------------------
echo "building the accession to symbol lookup ..."
gzip -dc "$DATASET" \
  | awk -F'\t' 'NR>1 && $11!="" && !seen[$11 FS $5]++ {print $11 FS $5}' \
  | awk -F'\t' '{
        if ($2 ~ /^ENSG/) { e[$1]=$2 }
        else { s[$1] = ($1 in s) ? s[$1] "|" $2 : $2 }
    }
    END { for (a in s) print a "\t" s[a]
          for (a in e) if (!(a in s)) print a "\t" e[a] }' \
  > "$WORK/uniprot_to_symbol.tsv"

# ---------------------------------------------------------------------------
# 5. Assemble. Every in-domain variant appears; those without homologous
#    pathogenic variation carry zero counts and empty accession lists.
# ---------------------------------------------------------------------------
echo "writing $FINAL ..."
{
    printf 'chrom\tpos\tclinvar_id\tref\talt\tclnsig\tclnrevstat\tclnhgvs\tgene_symbol\tuniprot_ac\tuniprot_pos\tpfam_id\tconsensus_pos\thomolog_missense_P_count\thomolog_missense_LP_count\tclinvar_P_accessions\tclinvar_LP_accessions\thomolog_benign_count\tsame_residue_PLP_count\tmetadome_url\n'
    awk -F'\t' -v OFS='\t' -v ANYFILE="$WORK/any_homologous.txt" '
      FILENAME==ARGV[1] { sym[$1]=$2; next }
      FILENAME==ARGV[2] { k=$4 SUBSEP $19 SUBSEP $20 SUBSEP $21 SUBSEP $22
                          P[k]=$23; LP[k]=$24; PA[k]=$25; LPA[k]=$26; next }
      # Benign variants per meta-position, and per residue within it, each
      # counted once per distinct ClinVar record.
      FILENAME==ARGV[3] { mk = $18 SUBSEP $19
                          rk = mk SUBSEP $16 SUBSEP $17
                          if (!((mk SUBSEP $4) in seen_pos)) { seen_pos[mk SUBSEP $4]; total[mk]++ }
                          if (!((rk SUBSEP $4) in seen_res)) { seen_res[rk SUBSEP $4]; here[rk]++ }
                          next }
      # Pathogenic variants per residue. No self-exclusion is needed: the
      # query variant is of uncertain significance and so cannot appear here.
      FILENAME==ARGV[4] { rk = $16 SUBSEP $17
                          if (!((rk SUBSEP $4) in seen_sp)) { seen_sp[rk SUBSEP $4]; sp[rk]++ }
                          next }
      { k = $4 SUBSEP $16 SUBSEP $17 SUBSEP $18 SUBSEP $19
        mk = $18 SUBSEP $19
        rk = mk SUBSEP $16 SUBSEP $17
        # Benign observations at this meta-position other than on this residue,
        # matching the self-exclusion the release applies to pathogenic counts.
        benign = (mk in total ? total[mk] : 0) - (rk in here ? here[rk] : 0)
        if (benign < 0) benign = 0
        # Any ClinVar missense evidence at this meta-position, benign included.
        if ((k in P) || benign > 0) any_ev[$4]
        if (!(k in P)) next
        same_residue = (($16 SUBSEP $17) in sp) ? sp[$16 SUBSEP $17] : 0
        print $1, $3, $4, $5, $6, $7, $8, $9,
              ($16 in sym ? sym[$16] : ""), $16, $17, $18, $19,
              P[k], LP[k], PA[k], LPA[k], benign, same_residue, $20 }
      END { n = 0; for (v in any_ev) n++; printf "%d\n", n > ANYFILE }
  ' "$WORK/uniprot_to_symbol.tsv" "$WORK/with_evidence.tsv" \
    "$WORK/benign_in_domain.tsv" "$WORK/same_protein_in_domain.tsv" \
    "$WORK/in_domain.tsv" \
    | sort -k1,1 -k2,2n
} > "$FINAL"

# ---------------------------------------------------------------------------
# 6. Summary, and the check that the coordinate conversion held.
#
# The release excludes a codon's own ClinVar records from its homologous
# counts, so a variant's own identifier must never appear in its own accession
# lists. A non-zero count means the coordinate shift in step 2 is wrong.
# ---------------------------------------------------------------------------
SELF=$(awk -F'\t' 'NR>1 && ($16 ~ $3 || $17 ~ $3)' "$FINAL" | wc -l | tr -d ' ')

echo
echo "  ${CLNSIG} missense SNVs    : $(wc -l < "$WORK/vus.bed" | tr -d ' ')"
echo "  rows in a Pfam domain       : $(wc -l < "$WORK/in_domain.tsv" | tr -d ' ')"
echo "  distinct VUS in a domain    : $(cut -f4 "$WORK/in_domain.tsv" | sort -u | wc -l | tr -d ' ')"
echo "  rows written (with P/LP)    : $(($(wc -l < "$FINAL" | tr -d ' ') - 1))"
echo "  distinct VUS written        : $(awk -F'\t' 'NR>1 {print $3}' "$FINAL" | sort -u | wc -l | tr -d ' ')"
echo "  benign missense SNVs        : $(wc -l < "$WORK/benign.bed" | tr -d ' ')"
echo "  pathogenic missense SNVs    : $(wc -l < "$WORK/pathogenic.bed" | tr -d ' ')"
echo "  distinct VUS, no same-residue: $(awk -F'\t' 'NR>1 && $19==0 {print $3}' "$FINAL" | sort -u | wc -l | tr -d ' ')"
echo "  distinct VUS, MetaDome-only : $(awk -F'\t' 'NR>1 && $18==0 && $19==0 {print $3}' "$FINAL" | sort -u | wc -l | tr -d ' ')"
echo "  distinct VUS, unconflicted  : $(awk -F'\t' 'NR>1 && $18==0 {print $3}' "$FINAL" | sort -u | wc -l | tr -d ' ')"
echo "  self-contamination          : $SELF (must be 0)"


if [ "$SELF" -ne 0 ]; then
    echo "WARNING: a variant appears among its own homologous accessions." >&2
fi

{
  printf 'genome_build\t%s\n' "$GENOME_BUILD"
  printf 'clinvar\t%s\n' "$(basename "$CLINVAR_VCF")"
  printf 'clnsig\t%s\n' "$CLNSIG"
  printf 'vus_missense_snvs\t%s\n' "$(wc -l < "$WORK/vus.bed" | tr -d ' ')"
  printf 'distinct_vus_in_domain\t%s\n' "$(cut -f4 "$WORK/in_domain.tsv" | sort -u | wc -l | tr -d ' ')"
  printf 'distinct_vus_any_homologous_missense\t%s\n' "$(cat "$WORK/any_homologous.txt")"
  printf 'distinct_vus_with_evidence\t%s\n' "$(awk -F'\t' 'NR>1 {print $3}' "$FINAL" | sort -u | wc -l | tr -d ' ')"
  printf 'benign_missense_snvs\t%s\n' "$(wc -l < "$WORK/benign.bed" | tr -d ' ')"
  printf 'distinct_vus_with_evidence_unconflicted\t%s\n' "$(awk -F'\t' 'NR>1 && $18==0 {print $3}' "$FINAL" | sort -u | wc -l | tr -d ' ')"
  printf 'pathogenic_missense_snvs\t%s\n' "$(wc -l < "$WORK/pathogenic.bed" | tr -d ' ')"
  printf 'distinct_vus_no_same_residue_evidence\t%s\n' "$(awk -F'\t' 'NR>1 && $19==0 {print $3}' "$FINAL" | sort -u | wc -l | tr -d ' ')"
  printf 'distinct_vus_metadome_only\t%s\n' "$(awk -F'\t' 'NR>1 && $18==0 && $19==0 {print $3}' "$FINAL" | sort -u | wc -l | tr -d ' ')"
} > "$OUT_DIR/vus_metadomain_${GENOME_BUILD}_summary.tsv"

if [ "$KEEP_INTERMEDIATES" -eq 0 ]; then
    rm -rf "$WORK"
else
    echo "  intermediates              : $WORK"
fi