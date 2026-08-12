#!/usr/bin/env python3
"""Plot the Figure 2 PWM similarity network.

Run from this directory:
    python E_plot_pwm_similarity_network.py
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
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap


SCRIPT_DIR = Path(__file__).resolve().parent
SSTAT_PATH = SCRIPT_DIR / "E_ssDNA-SELEX_sstat.txt"
OUTPUT_PREFIX = SCRIPT_DIR / "E_PWM_similarity_network"
EXCLUDED_PROTEIN_IDS = {"YP_001111258.1"}


def clean_label(x):
    """Clean PWM labels for plotting."""
    x = os.path.basename(str(x))
    x = x.replace(".txt", "")
    x = x.replace("_autoseed_pwm", "")
    x = x.replace("_pwm", "")
    x = x.replace("_PWM", "")
    parts = x.split("_")
    if parts and parts[-1] in {"A", "B"}:
        parts = parts[:-1]
    return "\n".join(parts)


def protein_prefix(x):
    return str(x).split("_PWM")[0]


def read_edges():
    sstat_similarity_df = pd.read_csv(SSTAT_PATH, sep="\t")

    required_columns = [
        "Protein_ID_A",
        "Protein_ID_B",
        "Similarity score (Ssum)",
    ]
    missing = [
        column for column in required_columns
        if column not in sstat_similarity_df.columns
    ]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    edge = sstat_similarity_df[
        ["Protein_ID_A", "Protein_ID_B", "Similarity score (Ssum)"]
    ].copy()
    edge = edge.rename(
        columns={
            "Protein_ID_A": "matrixA",
            "Protein_ID_B": "matrixB",
            "Similarity score (Ssum)": "Ssum",
        }
    )

    edge["Ssum"] = pd.to_numeric(edge["Ssum"], errors="coerce")
    edge = edge.dropna(subset=["matrixA", "matrixB", "Ssum"])
    edge["proteinA"] = edge["matrixA"].map(protein_prefix)
    edge["proteinB"] = edge["matrixB"].map(protein_prefix)
    edge = edge[
        (edge["proteinA"] == edge["proteinB"])
        & (~edge["proteinA"].isin(EXCLUDED_PROTEIN_IDS))
    ].copy()
    edge = edge.drop(columns=["proteinA", "proteinB"])

    if edge.empty:
        raise ValueError(f"No valid edges found in {SSTAT_PATH}")

    edge["matrixA"] = edge["matrixA"].map(clean_label)
    edge["matrixB"] = edge["matrixB"].map(clean_label)
    edge = edge[edge["matrixA"] != edge["matrixB"]].copy()

    edge[["matrixA", "matrixB"]] = pd.DataFrame(
        edge.apply(lambda row: sorted([row["matrixA"], row["matrixB"]]), axis=1).tolist(),
        index=edge.index,
    )
    edge = edge.groupby(["matrixA", "matrixB"], as_index=False)["Ssum"].max()

    return edge


def build_network(edge):
    G = nx.Graph()
    for _, row in edge.iterrows():
        G.add_edge(row["matrixA"], row["matrixB"], weight=row["Ssum"])
    return G


def component_layout(G):
    """Place disconnected protein components on a regular grid."""
    components = [
        sorted(component)
        for component in nx.connected_components(G)
    ]
    components = sorted(components, key=lambda component: component[0])

    positions = {}
    columns = 2
    x_gap = 6.0
    y_gap = 4.2
    scale = 1.4

    for index, component in enumerate(components):
        row = index // columns
        col = index % columns
        center_x = col * x_gap
        center_y = -row * y_gap
        subgraph = G.subgraph(component)
        local_pos = nx.circular_layout(subgraph, scale=scale)

        for node, (x_value, y_value) in local_pos.items():
            positions[node] = (x_value + center_x, y_value + center_y)

    return positions


def plot_network(G):
    pos = component_layout(G)

    fig, ax = plt.subplots(figsize=(18, 14))

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
        width=5,
        ax=ax,
    )

    nx.draw_networkx_nodes(
        G,
        pos,
        node_size=2800,
        node_color="#D9EAF7",
        edgecolors="#2B5C84",
        linewidths=2,
        ax=ax,
    )

    nx.draw_networkx_labels(
        G,
        pos,
        font_size=16,
        font_weight="bold",
        ax=ax,
    )

    sm = plt.cm.ScalarMappable(
        cmap=edge_cmap,
        norm=plt.Normalize(vmin=edge_vmin, vmax=edge_vmax),
    )
    sm.set_array([])

    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label("Similarity Score", fontsize=20, fontweight="bold")
    cbar.ax.tick_params(labelsize=16)
    for label in cbar.ax.get_yticklabels():
        label.set_fontweight("bold")

    ax.axis("off")
    ax.set_title(
        "PWM Similarity Network",
        fontsize=30,
        pad=18,
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
