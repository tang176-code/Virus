#!/usr/bin/env python3
"""Plot the HT-SELEX PWM similarity network.

Run from this directory:
    python C_plot_pwm_similarity_network.py
"""

import os
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Liberation Sans"]

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap


SCRIPT_DIR = Path(__file__).resolve().parent
SSTAT_PATH = SCRIPT_DIR / "C_HT-SELEX_sstat.txt"
OUTPUT_PREFIX = SCRIPT_DIR / "C_HT-SELEX_PWM_similarity_network"


def clean_label(x):
    """Clean PWM labels for plotting."""
    x = os.path.basename(str(x))
    x = x.replace(".txt", "")
    x = x.replace("_autoseed_pwm", "")
    x = x.replace("_pwm", "")
    x = x.replace("_PWM", "")
    x = "\n".join(x.split("_"))
    return x


def read_edges():
    sstat_similarity_df = pd.read_csv(SSTAT_PATH, sep="\t")

    required_columns = [
        "Protein_ID_A",
        "Protein_ID_B",
        "Similarity score",
    ]
    missing = [
        column for column in required_columns
        if column not in sstat_similarity_df.columns
    ]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    edge = sstat_similarity_df[
        ["Protein_ID_A", "Protein_ID_B", "Similarity score"]
    ].copy()
    edge = edge.rename(
        columns={
            "Protein_ID_A": "matrixA",
            "Protein_ID_B": "matrixB",
            "Similarity score": "Ssum",
        }
    )

    edge["Ssum"] = pd.to_numeric(edge["Ssum"], errors="coerce")
    edge = edge.dropna(subset=["matrixA", "matrixB", "Ssum"])

    if edge.empty:
        raise ValueError(f"No valid edges found in {SSTAT_PATH}")

    return edge


def build_network(edge):
    G = nx.Graph()
    for _, row in edge.iterrows():
        exp_a = clean_label(row["matrixA"])
        exp_b = clean_label(row["matrixB"])
        G.add_edge(exp_a, exp_b, weight=row["Ssum"])
    return G


def plot_network(G):
    pos = nx.spring_layout(G, k=5, seed=42, iterations=1500)

    fig, ax = plt.subplots(figsize=(32, 22))

    edge_colors = [G[u][v]["weight"] for u, v in G.edges()]
    edge_vmin = min(edge_colors)
    edge_vmax = max(edge_colors)

    if edge_vmin == edge_vmax:
        edge_vmin = edge_vmin - 0.01
        edge_vmax = edge_vmax + 0.01

    edge_cmap = LinearSegmentedColormap.from_list(
        "custom",
        ["#FFF3E2", "#E4003A"],
    )

    nx.draw_networkx_edges(
        G,
        pos,
        edge_color=edge_colors,
        edge_cmap=edge_cmap,
        edge_vmin=edge_vmin,
        edge_vmax=edge_vmax,
        width=6,
        ax=ax,
    )

    nx.draw_networkx_nodes(
        G,
        pos,
        node_size=3000,
        node_color="#D9EAF7",
        edgecolors="#2B5C84",
        linewidths=2,
        ax=ax,
    )

    nx.draw_networkx_labels(
        G,
        pos,
        font_size=20,
        font_weight="bold",
        ax=ax,
    )

    sm = plt.cm.ScalarMappable(
        cmap=edge_cmap,
        norm=plt.Normalize(vmin=edge_vmin, vmax=edge_vmax),
    )
    sm.set_array([])

    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label("Similarity Score", fontsize=24, fontweight="bold")
    cbar.ax.tick_params(labelsize=20)
    for label in cbar.ax.get_yticklabels():
        label.set_fontweight("bold")

    ax.axis("off")
    ax.set_title(
        "PWM Similarity Network",
        fontsize=38,
        pad=20,
        fontweight="bold",
    )

    plt.tight_layout()
    plt.savefig(f"{OUTPUT_PREFIX}.svg", format="svg", dpi=300, bbox_inches="tight")
    plt.savefig(f"{OUTPUT_PREFIX}.pdf", format="pdf", dpi=300, bbox_inches="tight")
    plt.savefig(f"{OUTPUT_PREFIX}.png", dpi=300, bbox_inches="tight")


def main():
    edge = read_edges()
    G = build_network(edge)
    plot_network(G)
    print(f"Read {len(edge)} edges from {SSTAT_PATH}")
    print(f"Network nodes: {G.number_of_nodes()}")
    print(f"Network edges: {G.number_of_edges()}")
    print(f"Saved {OUTPUT_PREFIX}.svg")
    print(f"Saved {OUTPUT_PREFIX}.pdf")
    print(f"Saved {OUTPUT_PREFIX}.png")


if __name__ == "__main__":
    main()
