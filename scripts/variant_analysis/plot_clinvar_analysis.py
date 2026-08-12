"""Figures for the ClinVar meta-domain analysis.

  Panel A follows the variants from ClinVar to the candidate set. Panel B shows
  how much homologous evidence each candidate carries, split by whether that
  evidence is pathogenic, likely pathogenic, or both.

  Reads the two files written by clinvar_analysis.sh for one genome build.
  """

import argparse
import os

import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

COLOUR_P_ONLY = "#a50f15"  # darkest: all homologues are Pathogenic
COLOUR_BOTH = "#d62728"  # MetaDome's pathogenic red
COLOUR_LP_ONLY = "#ff7f0e"  # MetaDome's likely-pathogenic orange

PATHOGENIC_ORDER = ["0", "1", "2", "3", "4", "5-9", "10+"]
TYPE_ORDER = ["Pathogenic only",
              "Pathogenic and likely pathogenic",
              "Likely pathogenic only"]
TYPE_COLOURS = {
    "Pathogenic only": COLOUR_P_ONLY,
    "Pathogenic and likely pathogenic": COLOUR_BOTH,
    "Likely pathogenic only": COLOUR_LP_ONLY,
}

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
    placement. The placement with the most Pathogenic homologues is kept,
    ties broken by likely pathogenic, so evidence shared between isoforms of the same protein is not counted twice.
    """
    frame = pd.read_csv(
        table, sep="\t", dtype={"clinvar_id": str},
        usecols=["clinvar_id", "homolog_missense_P_count", "homolog_missense_LP_count"],
    )
    frame = frame.sort_values(
        ["clinvar_id", "homolog_missense_P_count", "homolog_missense_LP_count"],
        ascending=[True, False, False])
    return frame.drop_duplicates("clinvar_id").reset_index(drop=True)


def pathogenic_bin(value):
    if value >= 10:
        return "10+"
    if value >= 5:
        return "5-9"
    return str(int(value))


def evidence_type(row):
    if row["homolog_missense_P_count"] > 0 and row["homolog_missense_LP_count"] > 0:
        return "Pathogenic and likely pathogenic"
    if row["homolog_missense_P_count"] > 0:
        return "Pathogenic only"
    return "Likely pathogenic only"


def panel_a(axis, summary, candidates):
    stages = ["Missense VUS\nin ClinVar",
              "Within a\nPfam domain",
              "With homologous\nP/LP missense"]
    values = [int(summary["vus_missense_snvs"]),
              int(summary["distinct_vus_in_domain"]),
              candidates]

    axis.bar(stages, [v / 1e6 for v in values],
             color=["#b0b0b0", "#7f7f7f", "#4d4d4d"], width=0.6)
    for index, value in enumerate(values):
        share = "" if index == 0 else "\n({:.1f}%)".format(100.0 * value / values[index - 1])
        axis.text(index, value / 1e6, "{:,}{}".format(value, share),
                  ha="center", va="bottom", fontsize=9)

    axis.set_ylabel("Variants (millions)")
    axis.set_ylim(0, values[0] / 1e6 * 1.25)
    axis.set_title("A", loc="left", fontweight="bold")
    sns.despine(ax=axis)


def panel_b(axis, grouped):
    grouped = grouped.copy()
    grouped["bin"] = grouped["homolog_missense_P_count"].apply(pathogenic_bin)
    grouped["type"] = grouped.apply(evidence_type, axis=1)

    counts = (grouped.groupby(["bin", "type"]).size()
              .unstack(fill_value=0)
              .reindex(index=PATHOGENIC_ORDER, fill_value=0)
              .reindex(columns=TYPE_ORDER, fill_value=0))

    total = len(grouped)
    share = counts / total * 100.0
    positions = list(range(len(PATHOGENIC_ORDER)))
    width = 0.26
    for offset, label in zip((-width, 0.0, width), TYPE_ORDER):
        axis.bar([p + offset for p in positions], share[label], width=width,
                 color=TYPE_COLOURS[label], label=label)
        for position, value in zip(positions, share[label]):
            if value == 0:
                axis.text(position + offset, 0.3, "0", ha="center", fontsize=7, color="#777777")

    axis.set_xticks(positions)
    axis.set_xticklabels(PATHOGENIC_ORDER)

    axis.axvline(2.5, color="#333333", linestyle="--", linewidth=1)
    axis.text(2.6, share.values.max() * 0.92,
              "\u2265 3 pathogenic observations:\nPM5 at moderate weight",
              fontsize=8, va="top", ha="left")
    axis.set_xlabel("Homologous pathogenic missense variants "
                    "at the equivalent meta-domain position")
    axis.set_ylabel("Share of candidate variants (%, n = {:,})".format(total))
    axis.legend(frameon=False, fontsize=9)
    axis.set_title("B", loc="left", fontweight="bold")
    sns.despine(ax=axis)

    return counts

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
    panel_a(axes[0], summary, len(grouped))
    counts = panel_b(axes[1], grouped)
    figure.suptitle("Meta-domain evidence for ClinVar variants of uncertain "
                    "significance ({})".format(args.genome_build), y=1.02)
    figure.tight_layout()

    stem = os.path.join(args.out_dir, "vus_metadomain_{}".format(args.genome_build))
    for extension in ("pdf", "png"):
        figure.savefig("{}.{}".format(stem, extension), dpi=300, bbox_inches="tight")
        print("wrote {}.{}".format(stem, extension))

    if args.counts:
        counts.to_csv(stem + "_panelB_counts.tsv", sep="\t")
        print("wrote {}_panelB_counts.tsv".format(stem))

    print("\n  candidates            : {:,}".format(len(grouped)))
    for threshold in (1, 2, 3, 5, 10):
        print("  >= {:<2} homologous Pathogenic : {:,}".format(
            threshold, int((grouped["homolog_missense_P_count"] >= threshold).sum())))


if __name__ == "__main__":
    main()