#!/usr/bin/env python3
"""Plot binding distributions across viral genome-size groups.

Run from this directory:
    python F-H_plot_binding_distribution.py
"""

from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Liberation Sans"]

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_PATH = SCRIPT_DIR / "F-H_BInding_distribution.csv"
OUTPUT_PREFIX = SCRIPT_DIR / "F-H_Binding_distribution_three_panel"
SOURCE_DATA_PATH = SCRIPT_DIR / "F-H_Binding_distribution_source_data.csv"

SIMILAR_LIST = [
    "YP_717939.1",
    "YP_010782993.1",
    "NP_043126.1",
    "YP_009111420.1",
    "YP_009553535.1",
    "AVM80381.1",
    "YP_008431134.1",
    "YP_717937.1",
    "YP_001111256.1",
]


def split_order(order):
    split_pid1 = "NP_957885.1"
    split_pid2 = "YP_005271187.1"

    if split_pid1 in order and split_pid2 in order:
        idx1 = order.index(split_pid1)
        idx2 = order.index(split_pid2)
        if idx1 < idx2:
            return [
                order[: idx1 + 1],
                order[idx1 + 1 : idx2 + 1],
                order[idx2 + 1 :],
            ]

    return [order, [], []]


def compact_kb_label(value):
    text = f"{value:.1f}"
    return text.rstrip("0").rstrip(".")


def write_source_data(data_df, ordersub_list, subtitles, protein_distance):
    source_rows = []

    for panel_index, (ordersub, subtitle) in enumerate(zip(ordersub_list, subtitles), start=1):
        if not ordersub:
            continue

        sub_df = data_df.set_index("Protein_ID").loc[ordersub].reset_index()
        x_pos = [i * protein_distance for i in range(len(ordersub))]

        for plot_order, (x_value, row) in enumerate(zip(x_pos, sub_df.to_dict("records")), start=1):
            genome_length_bp = row["Genome_length (bp)"]
            source_rows.append(
                {
                    "panel": subtitle.strip("()"),
                    "panel_index": panel_index,
                    "plot_order": plot_order,
                    "x_position": x_value,
                    "Protein_ID": row["Protein_ID"],
                    "Family": row["Family"],
                    "NBP_group": row["NBP_group"],
                    "highlight_in_similar_list": str(row["Protein_ID"]) in SIMILAR_LIST,
                    "TCPK_viral_genome": row["TCPK on Viral genome"],
                    "TCPK_host_chromosome": row["TCPK on host chromosome"],
                    "TCPK_host_mitochondrial_genome": row[
                        "TCPK on host Mitochondrial genome"
                    ],
                    "viral_genome_length_bp": genome_length_bp,
                    "viral_genome_size_kb": genome_length_bp / 1000.0,
                    "mitochondrial_genome_length_bp": row[
                        "Mitochondrial_genome_length(bp)"
                    ],
                    "FC_viral_vs_host_genome": row["FC_viral_vs_host_genome"],
                    "FC_MT_vs_viral": row["FC_MT_vs_viral"],
                }
            )

    pd.DataFrame(source_rows).to_csv(SOURCE_DATA_PATH, index=False)


def main():
    data_df = pd.read_csv(INPUT_PATH)
    numeric_columns = [
        "TCPK on Viral genome",
        "TCPK on host chromosome",
        "TCPK on host Mitochondrial genome",
        "Genome_length (bp)",
    ]
    for column in numeric_columns:
        data_df[column] = pd.to_numeric(data_df[column], errors="coerce")

    sorted_df = data_df[["Protein_ID", "Genome_length (bp)"]].sort_values(
        by=["Genome_length (bp)"],
        ascending=False,
    )
    order = sorted_df["Protein_ID"].dropna().astype(str).unique().tolist()
    ordersub_list = split_order(order)
    subtitles = ["(F)", "(G)", "(H)"]

    all_genome_kb = data_df["Genome_length (bp)"].dropna() / 1000.0
    right_ymax = all_genome_kb.max() * 1.10 if len(all_genome_kb) else None
    right_yticks = (
        list(range(0, int(right_ymax) + 26, 25))
        if right_ymax is not None
        else None
    )

    protein_distance = 1.6
    write_source_data(data_df, ordersub_list, subtitles, protein_distance)

    fig, axes = plt.subplots(1, 3, figsize=(30, 6), sharey=True)

    for ax, ordersub, subtitle in zip(axes, ordersub_list, subtitles):
        if not ordersub:
            ax.axis("off")
            continue

        sub_df = (
            data_df.set_index("Protein_ID")
            .loc[ordersub]
            .reset_index()
        )
        x_pos = [i * protein_distance for i in range(len(ordersub))]

        host_y = sub_df["TCPK on host chromosome"].tolist()
        viral_y = sub_df["TCPK on Viral genome"].tolist()
        mt_y = sub_df["TCPK on host Mitochondrial genome"].tolist()
        genome_kb = (sub_df["Genome_length (bp)"] / 1000.0).tolist()

        ax.plot(
            x_pos,
            host_y,
            color="black",
            linewidth=1.8,
            marker="o",
            markersize=6,
            label="Mean host genome binding frequency (/kb)",
            zorder=4,
        )
        ax.plot(
            x_pos,
            viral_y,
            color="red",
            linewidth=1.8,
            marker="o",
            markersize=6,
            label="Viral binding frequency (/kb)",
            zorder=4,
        )
        ax.plot(
            x_pos,
            mt_y,
            color="green",
            linewidth=1.8,
            marker="o",
            markersize=6,
            label="Mean mitochondrial genome binding frequency (/kb)",
            zorder=4,
        )

        ax.set_ylabel("Target Counts Per Kilobase", fontsize=10)
        ax.tick_params(axis="y", labelsize=10)
        ax.set_ylim(0, 300)
        ax.set_yticks(range(0, 301, 50))
        ax.tick_params(axis="x", rotation=90, labelsize=10)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(ordersub, rotation=90, fontsize=10)
        for label, pid in zip(ax.get_xticklabels(), ordersub):
            label.set_color("red" if str(pid) in SIMILAR_LIST else "black")

        ax2 = ax.twinx()
        ax2.scatter(x_pos, genome_kb, marker="D", s=36, color="tab:blue", zorder=3)
        ax2.set_ylabel("Viral genome size (kb)", fontsize=10)
        ax2.tick_params(axis="y", labelsize=10)
        if right_ymax is not None:
            ax2.set_ylim(0, right_ymax)
        if right_yticks is not None:
            ax2.set_yticks(right_yticks)

        genome_values = [g for g in genome_kb if pd.notna(g)]
        if genome_values:
            size_range = (
                f"(Genome size: {compact_kb_label(min(genome_values))}-"
                f"{compact_kb_label(max(genome_values))} kb)"
            )
        else:
            size_range = ""

        ax.set_title(
            f"{subtitle} Proteins: {ordersub[0]} ~ {ordersub[-1]} {size_range}",
            fontsize=14,
        )
        ax.grid(False)
        ax2.grid(False)

    black_line = Line2D(
        [0],
        [0],
        color="black",
        marker="o",
        linewidth=1,
        label="Mean host genome binding frequency (/kb)",
    )
    red_line = Line2D(
        [0],
        [0],
        color="red",
        marker="o",
        linewidth=1,
        label="Viral binding frequency (/kb)",
    )
    green_line = Line2D(
        [0],
        [0],
        color="green",
        marker="o",
        linewidth=1,
        label="Mean mitochondrial genome binding frequency (/kb)",
    )
    blue_dot = Line2D(
        [0],
        [0],
        marker="D",
        linestyle="None",
        color="tab:blue",
        label="Viral genome size (kb)",
    )

    fig.legend(
        handles=[black_line, red_line, green_line, blue_dot],
        loc="center left",
        bbox_to_anchor=(0.91, 0.5),
        borderaxespad=0,
        fontsize=10,
    )
    fig.subplots_adjust(left=0.04, right=0.88, bottom=0.34, top=0.86, wspace=0.35)

    plt.savefig(f"{OUTPUT_PREFIX}.pdf", format="pdf", dpi=300, bbox_inches="tight")
    plt.savefig(f"{OUTPUT_PREFIX}.png", dpi=300, bbox_inches="tight")

    print(f"Read {len(data_df)} rows from {INPUT_PATH}")
    print(f"Saved {SOURCE_DATA_PATH}")
    print(f"Saved {OUTPUT_PREFIX}.pdf")
    print(f"Saved {OUTPUT_PREFIX}.png")


if __name__ == "__main__":
    main()
