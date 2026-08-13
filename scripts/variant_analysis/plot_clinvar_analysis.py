"""Figures for the ClinVar meta-domain analysis.

Panel A follows missense variants of uncertain significance from ClinVar to
those with pathogenic variation at an equivalent meta-domain position. Panel B
divides those candidates into mutually exclusive evidence classes, ordered by
which evidence takes precedence under ACMG PM5.

Reads the table and summary written by clinvar_analysis.sh for one genome
build.
"""

import argparse
import os

import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

PM5_ORDER = ["Pathogenic (P/LP) also at the\nsame protein position (PM5 strong)",
             "Pathogenic (P/LP) at homologous\npositions only (PM5 moderate)",
             "Pathogenic and benign at\nhomologous positions (conflicting)"]
PM5_COLOURS = {PM5_ORDER[0]: "#d62728",  # MetaDome pathogenic red
               PM5_ORDER[1]: "#ff7f0e",  # MetaDome likely-pathogenic orange
               PM5_ORDER[2]: "#7f9bb5"}  # neutral: evidence in both directions


def read_summary(path):
    """The key/value counts that the candidate table cannot carry."""
    summary = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            key, _, value = line.rstrip("\n").partition("\t")
            summary[key] = value
    return summary


def per_variant(table):
    """Collapse placements to one row per ClinVar variant.

    A variant covered by several domain placements or isoforms appears once per
    placement. The placement carrying the strongest evidence is kept, so
    evidence shared between isoforms of one protein is not counted twice.
    """
    frame = pd.read_csv(
        table, sep="\t", dtype={"clinvar_id": str},
        usecols=["clinvar_id", "homolog_missense_P_count", "homolog_missense_LP_count",
                 "homolog_benign_count", "same_residue_PLP_count"])
    # Strongest placement per variant: same-residue evidence first, then
    # absence of conflicting benign evidence, then depth of pathogenic evidence.
    frame = frame.sort_values(
        ["clinvar_id", "same_residue_PLP_count", "homolog_benign_count",
         "homolog_missense_P_count"],
        ascending=[True, False, True, False])
    return frame.drop_duplicates("clinvar_id").reset_index(drop=True)


def pm5_class(row):
    """Evidence class, following the decision order of Gunning & Wright.

    Pathogenic variation at the variant's own residue is assessed first and
    outranks the paralogue route; a meta-position carrying benign as well as
    pathogenic observations is conflicting.
    """
    if row["same_residue_PLP_count"] > 0:
        return PM5_ORDER[0]
    if row["homolog_benign_count"] > 0:
        return PM5_ORDER[2]
    return PM5_ORDER[1]


def panel_a(axis, summary, grouped):
    stages = ["All",
        "Within a\nPfam domain",
        "With homologous\npathogenic (P/LP)\nmissense"]
    values = [int(summary["vus_missense_snvs"]),
        int(summary["distinct_vus_in_domain"]),
        len(grouped)]

    axis.bar(stages, [v / 1e6 for v in values],
        color=["#b0b0b0", "#7f7f7f", "#4d4d4d"], width=0.6)
    for index, value in enumerate(values):
        note = "" if index == 0 else "\n({:.1f}% of previous)".format(
            100.0 * value / values[index - 1])
        axis.text(index, value / 1e6, "{:,}{}".format(value, note),
                  ha="center", va="bottom", fontsize=8)

    axis.set_ylabel("Variants (millions)")
    axis.set_ylim(0, values[0] / 1e6 * 1.30)
    axis.set_title("From all missense VUS to candidates with\nhomologous evidence", loc="center", fontsize=10)
    axis.text(-0.14, 1.10, "A", transform=axis.transAxes,
              fontweight="bold", fontsize=12, va="top", ha="left")
    sns.despine(ax=axis)


def panel_b(axis, grouped):
    grouped = grouped.copy()
    grouped["class"] = grouped.apply(pm5_class, axis=1)
    counts = grouped["class"].value_counts().reindex(PM5_ORDER, fill_value=0)
    total = len(grouped)
    share = counts / total * 100.0

    positions = list(range(len(PM5_ORDER)))
    axis.bar(positions, share.values, width=0.6,
             color=[PM5_COLOURS[c] for c in PM5_ORDER])
    for position, value, count in zip(positions, share.values, counts.values):
        axis.text(position, value, "{:,}\n({:.1f}%)".format(count, value),
                  ha="center", va="bottom", fontsize=8)

    axis.set_ylim(0, share.values.max() * 1.20)
    axis.set_xticks(positions)
    axis.set_xticklabels(PM5_ORDER, fontsize=8)
    axis.set_ylabel("Share of variants (%)")
    axis.set_title("Evidence class of candidates (n = {:,})".format(total),
                   loc="center", fontsize=10)
    axis.text(-0.10, 1.10, "B", transform=axis.transAxes,
              fontweight="bold", fontsize=12, va="top", ha="left")
    sns.despine(ax=axis)

    return counts.to_frame("variants")

def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--table", required=True,
                        help="vus_metadomain_<build>.tsv")
    parser.add_argument("--summary", required=True,
                        help="vus_metadomain_<build>_summary.tsv")
    parser.add_argument("--genome-build", required=True)
    parser.add_argument("--out-dir", default=".")
    parser.add_argument("--counts", action="store_true",
                        help="also write the panel B counts as TSV")
    args = parser.parse_args()

    sns.set_theme(style="ticks", context="paper")

    summary = read_summary(args.summary)
    grouped = per_variant(args.table)

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2),
                                gridspec_kw={"width_ratios": [1, 1.6]})
    panel_a(axes[0], summary, grouped)
    counts = panel_b(axes[1], grouped)
    figure.suptitle("Meta-domain PM5 evidence for ClinVar missense variants of uncertain "
                    "significance ({})".format(args.genome_build), y=1.02)
    figure.tight_layout()

    stem = os.path.join(args.out_dir, "vus_metadomain_{}".format(args.genome_build))
    for extension in ("pdf", "png"):
        figure.savefig("{}.{}".format(stem, extension), dpi=300, bbox_inches="tight")
        print("wrote {}.{}".format(stem, extension))

    if args.counts:
        counts.to_csv(stem + "_panelB_counts.tsv", sep="\t")
        print("wrote {}_panelB_counts.tsv".format(stem))


if __name__ == "__main__":
    main()