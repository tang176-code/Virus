#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Summarize the mtDNA PWM base-label shuffle null analysis.

This script is the compact analysis code for the GitHub version. The upstream
FIMO step counts mtDNA motif hits for:

1. the original viral PWM; and
2. 1,000 matched shuffled PWMs.

For each shuffled PWM, A/C/G/T labels are randomly permuted within every PWM
column. This preserves motif width and per-position information content, while
disrupting the original base preference.

The script starts from those observed and shuffled hit counts, then calculates
fold enrichment, empirical P values, null intervals, z-scores, and BH-FDR.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path(__file__).resolve().parent
DEFAULT_OBSERVED = BASE / "mtdna_base_label_shuffle_null_observed_counts.csv"
DEFAULT_SHUFFLED = BASE / "mtdna_base_label_shuffle_null_shuffled_counts.csv"
DEFAULT_OUTPUT = BASE / "mtdna_base_label_shuffle_null_source_data.csv"

EXCLUDED_NON_HUMAN_INFECTING = {
    "AVM80381.1",
    "NP_671816.1",
    "YP_010782993.1",
}

OUTPUT_COLUMNS = [
    "Protein_ID",
    "genome_length",
    "Genome_length_group",
    "PWM_width",
    "FIMO_p_threshold",
    "observed_mtDNA_hits",
    "shuffled_null_mean",
    "shuffled_null_SD",
    "empirical_null_2.5_percentile",
    "empirical_null_97.5_percentile",
    "fold_enrichment",
    "z_score",
    "n_shuffled_ge_observed",
    "empirical_P_over",
    "n_shuffled_PWMs",
    "n_identity_shuffles_rejected",
    "BH_adjusted_empirical_P",
    "significant_BH_FDR_0.05",
]

COLUMN_DESCRIPTIONS = {
    "Protein_ID": "Database accession number of the viral nucleic acid-binding protein.",
    "genome_length": "Length of the corresponding viral genome, in base pairs.",
    "Genome_length_group": "Viral genome-size category: Large, Medium, or Small.",
    "PWM_width": "Width of the PWM, corresponding to the number of nucleotide positions in the motif.",
    "FIMO_p_threshold": "FIMO P-value filtering threshold used to retain motif hits.",
    "observed_mtDNA_hits": "Number of mtDNA motif hits detected using the original PWM at the indicated filtering threshold.",
    "shuffled_null_mean": "Mean mtDNA motif-hit count obtained from the matched base-label-shuffled PWMs.",
    "shuffled_null_SD": "Standard deviation of mtDNA motif-hit counts in the shuffled-PWM null distribution.",
    "empirical_null_2.5_percentile": "2.5th percentile of the shuffled-PWM null distribution of mtDNA motif-hit counts.",
    "empirical_null_97.5_percentile": "97.5th percentile of the shuffled-PWM null distribution of mtDNA motif-hit counts.",
    "fold_enrichment": "Observed mtDNA hit count divided by the shuffled-PWM null mean; values above 1 indicate enrichment relative to the null mean.",
    "z_score": "Standardized difference, calculated as (observed hits - null mean) / null SD.",
    "n_shuffled_ge_observed": "Number of shuffled PWMs with mtDNA hit counts greater than or equal to the observed count.",
    "empirical_P_over": "One-sided empirical P value, calculated as (1 + number of shuffled PWMs with hits >= observed hits) / (1 + total shuffled PWMs).",
    "n_shuffled_PWMs": "Number of matched shuffled PWMs generated for each original PWM.",
    "n_identity_shuffles_rejected": "Number of identity shuffle attempts that were rejected and regenerated.",
    "BH_adjusted_empirical_P": "BH-adjusted empirical P value, calculated separately across proteins at each FIMO threshold.",
    "significant_BH_FDR_0.05": "Indicates whether the result remained significant after BH correction at FDR < 0.05.",
}


def base_label_shuffle_pwm(pwm: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Return one base-label-shuffled PWM.

    Rows are PWM positions and columns are A/C/G/T probabilities. Each row keeps
    the same four probability values, but reassigns them to nucleotide labels.
    """
    shuffled = np.empty_like(pwm)
    for i, row in enumerate(pwm):
        shuffled[i] = row[rng.permutation(4)]
    return shuffled


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    """BH-adjust P values while preserving the input index."""
    p = p_values.astype(float).to_numpy()
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    adjusted_ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty(n, dtype=float)
    adjusted[order] = np.minimum(adjusted_ranked, 1.0)
    return pd.Series(adjusted, index=p_values.index)


def normalize_observed_table(path: Path) -> pd.DataFrame:
    """Read observed original-PWM counts and metadata."""
    observed = pd.read_csv(path)
    observed = observed.rename(
        columns={
            "width": "PWM_width",
            "p_threshold": "FIMO_p_threshold",
            "observed": "observed_mtDNA_hits",
            "n_identity_skipped": "n_identity_shuffles_rejected",
        }
    )

    keep = [
            "Protein_ID",
        "genome_length",
        "Genome_length_group",
        "PWM_width",
        "FIMO_p_threshold",
        "observed_mtDNA_hits",
        "n_identity_shuffles_rejected",
    ]
    missing = sorted(set(keep) - set(observed.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {missing}")

    observed = observed[keep].copy()
    observed = observed[~observed["Protein_ID"].isin(EXCLUDED_NON_HUMAN_INFECTING)]
    return observed


def summarize_null_counts(observed: pd.DataFrame, shuffled_path: Path) -> pd.DataFrame:
    """Calculate per-protein statistics from shuffled-PWM mtDNA hit counts."""
    shuffled = pd.read_csv(shuffled_path).rename(
        columns={
            "p_threshold": "FIMO_p_threshold",
            "null_count": "shuffled_mtDNA_hits",
        }
    )
    shuffled = shuffled[~shuffled["Protein_ID"].isin(EXCLUDED_NON_HUMAN_INFECTING)]

    keys = ["Protein_ID", "FIMO_p_threshold"]
    merged = shuffled.merge(
        observed[keys + ["observed_mtDNA_hits"]],
        on=keys,
        how="inner",
        validate="many_to_one",
    )
    merged["shuffled_ge_observed"] = (
        merged["shuffled_mtDNA_hits"] >= merged["observed_mtDNA_hits"]
    ).astype(int)

    null_stats = (
        merged.groupby(keys)["shuffled_mtDNA_hits"]
        .agg(
            shuffled_null_mean="mean",
            shuffled_null_SD=lambda x: x.std(ddof=1),
            **{
                "empirical_null_2.5_percentile": lambda x: np.percentile(x, 2.5),
                "empirical_null_97.5_percentile": lambda x: np.percentile(x, 97.5),
            },
            n_shuffled_PWMs="size",
        )
        .reset_index()
    )
    ge_counts = (
        merged.groupby(keys)["shuffled_ge_observed"]
        .sum()
        .rename("n_shuffled_ge_observed")
        .reset_index()
    )

    result = observed.merge(null_stats, on=keys, how="inner", validate="one_to_one")
    result = result.merge(ge_counts, on=keys, how="inner", validate="one_to_one")
    result["fold_enrichment"] = (
        result["observed_mtDNA_hits"] / result["shuffled_null_mean"].replace(0, np.nan)
    )
    result["z_score"] = (
        (result["observed_mtDNA_hits"] - result["shuffled_null_mean"])
        / result["shuffled_null_SD"].replace(0, np.nan)
    )
    result["empirical_P_over"] = (
        (1 + result["n_shuffled_ge_observed"]) / (1 + result["n_shuffled_PWMs"])
    )

    result["BH_adjusted_empirical_P"] = np.nan
    for threshold, idx in result.groupby("FIMO_p_threshold").groups.items():
        result.loc[idx, "BH_adjusted_empirical_P"] = benjamini_hochberg(
            result.loc[idx, "empirical_P_over"]
        )
    result["significant_BH_FDR_0.05"] = result["BH_adjusted_empirical_P"] <= 0.05

    result["Genome_length_group"] = (
        result["Genome_length_group"].astype(str).replace({"Mid": "Medium"})
    )
    group_order = pd.CategoricalDtype(["Large", "Medium", "Small"], ordered=True)
    result["Genome_length_group"] = result["Genome_length_group"].astype(group_order)
    result = result.sort_values(
        ["FIMO_p_threshold", "Genome_length_group", "genome_length"],
        ascending=[False, True, False],
    )
    result["Genome_length_group"] = result["Genome_length_group"].astype(str)
    return result[OUTPUT_COLUMNS]


def write_csv_with_metadata(data: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("# CSV column descriptions\n")
        handle.write(
            "# Lines beginning with # are metadata comments. "
            "To read the data table with pandas, use pd.read_csv(path, comment='#').\n"
        )
        for column in OUTPUT_COLUMNS:
            handle.write(f"# {column}: {COLUMN_DESCRIPTIONS[column]}\n")
        handle.write("#\n")
        data.to_csv(handle, index=False)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Recalculate mtDNA PWM base-label shuffle null summary statistics."
    )
    parser.add_argument("--observed", type=Path, default=DEFAULT_OBSERVED)
    parser.add_argument("--shuffled", type=Path, default=DEFAULT_SHUFFLED)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    observed = normalize_observed_table(args.observed)
    summary = summarize_null_counts(observed, args.shuffled)
    write_csv_with_metadata(summary, args.out)

    print(f"Wrote {args.out}")
    print(f"Proteins: {summary['Protein_ID'].nunique()}")
    print(f"Rows: {len(summary)}")
    print(f"Excluded proteins: {', '.join(sorted(EXCLUDED_NON_HUMAN_INFECTING))}")


if __name__ == "__main__":
    main()
