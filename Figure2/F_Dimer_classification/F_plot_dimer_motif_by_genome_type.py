#!/usr/bin/env python3
"""Build dimer source data and plot motif class counts by genome type."""

from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Liberation Sans"]

import matplotlib.pyplot as plt
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
CLASSIFICATION_PATH = SCRIPT_DIR / "classification.tsv"
VIRAL_INFO_PATH = SCRIPT_DIR / "viral_information.csv"
SOURCE_DATA_PATH = SCRIPT_DIR / "F_dimer_motif_by_genome_type_source_data.csv"
COUNT_DATA_PATH = SCRIPT_DIR / "F_dimer_motif_by_genome_type_counts.csv"
OUTPUT_PREFIX = SCRIPT_DIR / "F_dimer_motif_by_genome_type_stacked_bar"

GENOME_TYPE_LABELS = {
    "ds-DNA": "dsDNA",
    "ss-DNA": "ssDNA",
}
GENOME_TYPE_ORDER = ["dsDNA", "ssDNA"]
MOTIF_ORDER = ["Monomer", "Dimer"]
MOTIF_COLORS = {
    "Monomer": "#B9C0C7",
    "Dimer": "#2FA4D7",
}


def protein_id_from_filename(filename):
    value = str(filename)
    if value.startswith("HT-SELEX_") and value.endswith("_pwm.txt"):
        return value.removeprefix("HT-SELEX_").removesuffix("_pwm.txt")
    return pd.NA


def build_source_data():
    classification = pd.read_csv(CLASSIFICATION_PATH, sep="\t")
    viral_info = pd.read_csv(VIRAL_INFO_PATH)

    classification["Protein_ID"] = classification["filename"].map(
        protein_id_from_filename
    )
    classification["Motif"] = classification["classification"].map(
        {"doublet": "Dimer", "singlet": "Monomer"}
    )

    required = {"filename", "Protein_ID", "Motif", "p_value", "q_bh"}
    missing = required.difference(classification.columns)
    if missing:
        raise ValueError(f"Missing required classification columns: {sorted(missing)}")
    if classification["Protein_ID"].isna().any():
        bad = classification.loc[classification["Protein_ID"].isna(), "filename"]
        raise ValueError(f"Could not parse Protein_ID from filenames: {bad.tolist()}")

    viral_info = viral_info.rename(
        columns={
            "Protein ID": "Protein_ID",
            "Genome Type": "Genome_Type",
            "Genome_length(bp)": "Genome_length_bp",
        }
    )
    metadata_columns = [
        "Protein_ID",
        "Genome_Type",
        "Accession",
        "Organism_Name",
        "Family",
        "Genome_length_bp",
    ]

    merged = classification.merge(
        viral_info[metadata_columns],
        on="Protein_ID",
        how="left",
        validate="one_to_one",
    )
    if merged["Genome_Type"].isna().any():
        missing_ids = merged.loc[merged["Genome_Type"].isna(), "Protein_ID"].tolist()
        raise ValueError(f"Missing Genome_Type for Protein_ID: {missing_ids}")

    merged["Genome_Type"] = merged["Genome_Type"].map(GENOME_TYPE_LABELS)
    merged = merged.loc[merged["Genome_Type"].isin(GENOME_TYPE_ORDER)].copy()
    merged["Genome_Type"] = pd.Categorical(
        merged["Genome_Type"],
        categories=GENOME_TYPE_ORDER,
        ordered=True,
    )
    merged["Motif"] = pd.Categorical(
        merged["Motif"],
        categories=MOTIF_ORDER,
        ordered=True,
    )
    merged = merged.sort_values(["Genome_Type", "Motif", "Protein_ID"]).reset_index(
        drop=True
    )
    merged.insert(0, "plot", range(1, len(merged) + 1))

    source_data = merged[
        [
            "plot",
            "Protein_ID",
            "Motif",
            "Genome_Type",
            "Accession",
            "Organism_Name",
            "Family",
            "Genome_length_bp",
            "p_value",
            "q_bh",
            "filename",
        ]
    ].copy()

    counts = (
        source_data.groupby(["Genome_Type", "Motif"], observed=False)
        .size()
        .unstack(fill_value=0)
        .reindex(index=GENOME_TYPE_ORDER, columns=MOTIF_ORDER, fill_value=0)
    )

    return source_data, counts


def plot_counts(counts):
    fig, ax = plt.subplots(figsize=(4.2, 3.2))

    x_values = range(len(counts.index))
    bottom = [0] * len(counts.index)
    for motif in MOTIF_ORDER:
        values = counts[motif].tolist()
        ax.bar(
            x_values,
            values,
            bottom=bottom,
            color=MOTIF_COLORS[motif],
            edgecolor="black",
            linewidth=0.4,
            label=motif,
            width=0.62,
        )
        for x_value, value, base in zip(x_values, values, bottom):
            if value > 0:
                ax.text(
                    x_value,
                    base + value / 2,
                    str(int(value)),
                    ha="center",
                    va="center",
                    fontsize=9,
                    fontweight="bold",
                )
        bottom = [base + value for base, value in zip(bottom, values)]

    ax.set_xticks(list(x_values))
    ax.set_xticklabels(counts.index.tolist(), fontsize=10)
    ax.set_ylabel("Number of proteins", fontsize=10)
    ax.set_xlabel("Genome type", fontsize=10)
    ax.tick_params(axis="y", labelsize=10)
    ax.set_ylim(0, max(bottom) + 2)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_linewidth(0.6)

    plt.tight_layout()
    plt.savefig(f"{OUTPUT_PREFIX}.pdf", format="pdf", dpi=300, bbox_inches="tight")
    plt.savefig(f"{OUTPUT_PREFIX}.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{OUTPUT_PREFIX}.svg", format="svg", dpi=300, bbox_inches="tight")


def main():
    source_data, counts = build_source_data()
    source_data.to_csv(SOURCE_DATA_PATH, index=False)
    counts.reset_index().rename(columns={"index": "Genome_Type"}).to_csv(
        COUNT_DATA_PATH,
        index=False,
    )
    plot_counts(counts)

    print(f"Saved {SOURCE_DATA_PATH}")
    print(f"Saved {COUNT_DATA_PATH}")
    print(f"Saved {OUTPUT_PREFIX}.pdf")
    print(f"Saved {OUTPUT_PREFIX}.png")
    print(f"Saved {OUTPUT_PREFIX}.svg")


if __name__ == "__main__":
    main()
