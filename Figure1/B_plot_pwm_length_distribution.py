#!/usr/bin/env python3
"""Plot the PWM length distribution for Figure 1.

Run from this directory:
    python B_plot_pwm_length_distribution.py
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
import seaborn as sns


SCRIPT_DIR = Path(__file__).resolve().parent
PWM_DIR = SCRIPT_DIR / "PWM"
LENGTH_TABLE = SCRIPT_DIR / "B_PWM_length.txt"
OUTPUT_PREFIX = SCRIPT_DIR / "B_PWM_length_distribution"


def protein_id_from_pwm_name(path):
    name = path.name
    prefix = "HT-SELEX_"
    suffix = "_pwm.txt"
    if not name.startswith(prefix) or not name.endswith(suffix):
        raise ValueError(f"Unexpected PWM file name: {name}")
    return name[len(prefix) : -len(suffix)]


def pwm_length(path):
    with path.open(encoding="utf-8-sig") as handle:
        header = handle.readline().strip().split("\t")
    if len(header) < 2 or header[0].lower() != "base":
        raise ValueError(f"Unexpected PWM header in {path.name}")
    return len(header) - 1


def read_pwm_lengths():
    records = []
    for pwm_file in sorted(PWM_DIR.glob("HT-SELEX_*_pwm.txt")):
        records.append(
            {
                "Protein_ID": protein_id_from_pwm_name(pwm_file),
                "motif_length": pwm_length(pwm_file),
                "PWM_file": pwm_file.name,
            }
        )
    if not records:
        raise FileNotFoundError(f"No PWM files found in {PWM_DIR}")
    return pd.DataFrame(records)


def verify_length_table(pwm_df):
    length_df = pd.read_csv(LENGTH_TABLE, sep="\t")
    merged = pwm_df.merge(
        length_df,
        on=["Protein_ID", "PWM_file"],
        how="outer",
        suffixes=("_from_pwm", "_from_table"),
        indicator=True,
    )

    unmatched = merged[merged["_merge"] != "both"]
    if not unmatched.empty:
        raise ValueError(
            "PWM files and PWM_length.txt do not contain the same records:\n"
            + unmatched.to_string(index=False)
        )

    mismatched = merged[
        merged["motif_length_from_pwm"] != merged["motif_length_from_table"]
    ]
    if not mismatched.empty:
        raise ValueError(
            "Motif lengths calculated from PWM do not match PWM_length.txt:\n"
            + mismatched.to_string(index=False)
        )


def plot_length_distribution(pwm_df):
    motif_lengths = pwm_df["motif_length"].astype(int)
    length_counts = motif_lengths.value_counts().sort_index()

    plt.figure(figsize=(6, 6))
    ax = sns.countplot(
        x=motif_lengths,
        order=length_counts.index.tolist(),
        color="#2FA4D7",
    )

    ax.set_xlabel("PWM length (bp)", fontsize=16, fontweight="bold")
    ax.set_ylabel("Number of models", fontsize=16, fontweight="bold")

    for p in ax.patches:
        height = int(p.get_height())
        ax.annotate(
            f"{height}",
            (p.get_x() + p.get_width() / 2.0, height),
            ha="center",
            va="center",
            xytext=(0, 5),
            textcoords="offset points",
            fontsize=12,
            fontweight="bold",
        )

    for label in ax.get_xticklabels():
        label.set_fontsize(14)
        label.set_fontweight("bold")
    for label in ax.get_yticklabels():
        label.set_fontsize(14)
        label.set_fontweight("bold")

    ax.set_facecolor("white")
    ax.grid(False)
    ax.set_axisbelow(False)
    sns.despine(ax=ax)

    plt.tight_layout()
    plt.savefig(f"{OUTPUT_PREFIX}.pdf", bbox_inches="tight")
    plt.savefig(f"{OUTPUT_PREFIX}.png", dpi=300, bbox_inches="tight")


def main():
    pwm_df = read_pwm_lengths()
    verify_length_table(pwm_df)
    plot_length_distribution(pwm_df)
    print(f"Read {len(pwm_df)} PWM files from {PWM_DIR}")
    print(f"Saved {OUTPUT_PREFIX}.pdf")
    print(f"Saved {OUTPUT_PREFIX}.png")


if __name__ == "__main__":
    main()
