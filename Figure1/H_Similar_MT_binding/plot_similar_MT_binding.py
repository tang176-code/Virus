#!/usr/bin/env python3
"""Plot mitochondrial and viral binding metrics for similar HT-SELEX PWMs.

Run from this directory:
    python plot_similar9_MT_binding.py
"""

from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Liberation Sans"]

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D


SCRIPT_DIR = Path(__file__).resolve().parent
MT_BINDING_PATH = SCRIPT_DIR / "similar7_MT_binding.csv"
SUMMARY_PATH = SCRIPT_DIR / "similar7_MT_binding_summary.csv"
OUTPUT_PREFIX = SCRIPT_DIR / "similar7virus_MT_dot_order_genomelength_zhexian"

SIMILAR_LIST = [
   yang
]


def prepare_mt_binding(mt_df):
    mt_df = mt_df[mt_df["Protein_ID"].isin(SIMILAR_LIST)].copy()

    fragment_pos = mt_df["fragment"].astype(str).str.split("REF_").str[-1]
    mt_df["start"] = fragment_pos.str.split("_").str[0].astype(int)
    mt_df["stop"] = fragment_pos.str.split("_").str[-1].astype(int)
    mt_df["fimo_scan_length"] = mt_df["stop"] - mt_df["start"] + 1

    mt_df["hit_count"] = pd.to_numeric(mt_df["hit_count"], errors="coerce")
    mt_df["fimo_scan_length"] = pd.to_numeric(
        mt_df["fimo_scan_length"],
        errors="coerce",
    )
    mt_df["Target_Counts_Per_Kilobase"] = (
        mt_df["hit_count"] / mt_df["fimo_scan_length"] * 1000
    )
    return mt_df


def main():
    mt_df = prepare_mt_binding(pd.read_csv(MT_BINDING_PATH))
    summary_df = pd.read_csv(SUMMARY_PATH)
    summary_df = summary_df[summary_df["Protein_ID"].isin(SIMILAR_LIST)].copy()

    summary_df = summary_df.sort_values("plot_order")
    sorted_list = summary_df["Protein_ID"].tolist()
    mean_by_protein = (
        mt_df.groupby("Protein_ID")["Target_Counts_Per_Kilobase"].mean()
    )
    summary_df["mean_host_genome_binding_frequency_per_kb"] = (
        summary_df["Protein_ID"].map(mean_by_protein)
    )

    fig_w, fig_h = 40 / 25.4, 25 / 25.4
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    sns.boxplot(
        x="Protein_ID",
        y="Target_Counts_Per_Kilobase",
        data=mt_df,
        order=sorted_list,
        color="#90C67C",
        showfliers=False,
        ax=ax,
        linewidth=0.25,
    )

    x_pos = list(range(len(sorted_list)))
    genome_kb_list = summary_df["genome_size_kb"].tolist()
    mean_list = summary_df["mean_host_genome_binding_frequency_per_kb"].tolist()
    vfreq_list = summary_df["viral_binding_frequency_per_kb"].tolist()

    line_x_mean = [x for x, y in zip(x_pos, mean_list) if pd.notna(y)]
    line_y_mean = [y for y in mean_list if pd.notna(y)]
    line_x_vfreq = [x for x, y in zip(x_pos, vfreq_list) if pd.notna(y)]
    line_y_vfreq = [y for y in vfreq_list if pd.notna(y)]

    ax.plot(
        line_x_vfreq,
        line_y_vfreq,
        color="red",
        linewidth=0.5,
        marker="o",
        markersize=1,
        label="Viral binding frequency (/kb)",
        zorder=1,
    )
    ax.plot(
        line_x_mean,
        line_y_mean,
        color="green",
        linewidth=0.5,
        marker="o",
        markersize=1,
        label="Mean host genome binding frequency (/kb)",
        zorder=0.3,
    )

    ax2 = ax.twinx()
    plot_x = [x for x, y in zip(x_pos, genome_kb_list) if pd.notna(y)]
    plot_y = [y for y in genome_kb_list if pd.notna(y)]
    ax2.scatter(plot_x, plot_y, marker="D", s=4, color="tab:blue", zorder=0.1)
    ax2.set_ylabel("Viral genome size (kb)", fontsize=8, labelpad=1)
    ax2.tick_params(axis="y", labelsize=8, width=0.5, length=0.5)

    ax.set_ylabel("Target Counts Per Kilobase", fontsize=8, labelpad=1)
    ax.tick_params(axis="y", labelsize=8, width=0.5, length=2.5)
    ax.set_ylim(0, 300)
    ax.set_yticks(range(0, 301, 50))
    ax.tick_params(axis="x", rotation=90, labelsize=8, width=0.5, length=2.5)

    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    for spine in ax2.spines.values():
        spine.set_linewidth(0.5)

    green_patch = mpatches.Patch(
        color="#90C67C",
        label="Mitochondrial genome binding (distribution)",
    )
    green_line = Line2D(
        [0],
        [0],
        color="green",
        marker="o",
        linewidth=0.5,
        markersize=3,
        label="Mean host genome binding frequency (/kb)",
    )
    red_line = Line2D(
        [0],
        [0],
        color="red",
        marker="o",
        linewidth=0.5,
        markersize=3,
        label="Viral binding freq (/kb)",
    )
    blue_dot = Line2D(
        [0],
        [0],
        marker="D",
        linestyle="None",
        color="tab:blue",
        label="Viral genome size (kb)",
        markersize=3,
    )

    ax.legend(
        handles=[green_patch, green_line, red_line, blue_dot],
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        borderaxespad=0,
        fontsize=8,
        frameon=False,
    )

    ax.grid(False)
    ax2.grid(False)
    plt.tight_layout(pad=0.12)

    plt.savefig(f"{OUTPUT_PREFIX}.pdf", format="pdf", dpi=300, bbox_inches="tight")
    plt.savefig(f"{OUTPUT_PREFIX}.png", dpi=300, bbox_inches="tight")

    print(f"Read {len(mt_df)} MT binding rows from {MT_BINDING_PATH}")
    print(f"Read {len(summary_df)} summary rows from {SUMMARY_PATH}")
    print(f"Saved {OUTPUT_PREFIX}.pdf")
    print(f"Saved {OUTPUT_PREFIX}.png")


if __name__ == "__main__":
    main()
