#!/usr/bin/env python3
"""Reproduce Similar_list three-context pooled hypergeometric enrichment outputs.

Inputs expected by default, relative to this script:
- data/Binding.csv
- data/merge_df_similar_list.csv.xz
- data/PWM_fimo_result/<Protein_ID>_fimo.tsv

Outputs written by default to results/ next to this script:
- Similar_list_three_context_hypergeometric_per_protein_results.csv
- Similar_list_three_context_hypergeometric_pooled_results.csv
- Similar_list_three_context_pooled_found_ids.csv
- Similar_list_three_context_missing_ids.csv
- Similar_list_three_context_hypergeometric_pooled_enrichment.pdf/.png
- Similar_list_three_context_hypergeometric_pooled_observed_expected_heatmap.pdf/.png
"""

from pathlib import Path
import argparse
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = SCRIPT_DIR / "results"
DEFAULT_BINDING_CSV = SCRIPT_DIR / "data" / "Binding.csv"
DEFAULT_MERGE_DF_CSV = SCRIPT_DIR / "data" / "merge_df_similar_list.csv.xz"
DEFAULT_FIMO_DIR = SCRIPT_DIR / "data" / "PWM_fimo_result"

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

CONTEXTS = ["Viral genome", "Host mitochondrial genome", "Host nuclear genome"]
CONTEXT_LABELS = {
    "Viral genome": "Viral genome",
    "Host mitochondrial genome": "mtDNA",
    "Host nuclear genome": "Nuclear genome",
}
COLORS = {
    "Viral genome": "#4C78A8",
    "Host mitochondrial genome": "#D55E00",
    "Host nuclear genome": "#54A24B",
}
ALPHA = 0.05


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Reproduce the Similar_list three-context pooled hypergeometric "
            "enrichment CSV, PDF, and PNG outputs from raw FIMO counts."
        )
    )
    parser.add_argument("--binding-csv", type=Path, default=DEFAULT_BINDING_CSV, help="Path to Binding.csv.")
    parser.add_argument("--merge-df", type=Path, default=DEFAULT_MERGE_DF_CSV, help="Path to merge_df.csv.")
    parser.add_argument("--fimo-dir", type=Path, default=DEFAULT_FIMO_DIR, help="Directory containing *_fimo.tsv files.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_RAW_DIR, help="Directory for generated CSV/PDF/PNG files.")
    return parser.parse_args()


def validate_inputs(binding_csv, merge_df_csv, fimo_dir):
    missing = []
    if not binding_csv.exists():
        missing.append(str(binding_csv))
    if not merge_df_csv.exists():
        missing.append(str(merge_df_csv))
    if not fimo_dir.exists():
        missing.append(str(fimo_dir))
    if missing:
        raise FileNotFoundError("Required input path(s) not found:\n" + "\n".join(missing))


def parse_scan_length(fragment_series):
    coords = fragment_series.str.split("REF_", n=1).str[-1].str.rsplit("_", n=1, expand=True)
    start = coords[0].astype(np.int64)
    stop = coords[1].astype(np.int64)
    return stop - start + 1


def summarize_merge_df(merge_df_csv, target_protein_ids):
    usecols = ["fragment", "hit_count", "chr", "virus_id", "Protein_ID"]
    dtype = {
        "fragment": "string",
        "hit_count": "int64",
        "chr": "string",
        "virus_id": "string",
        "Protein_ID": "string",
    }
    target_protein_ids = set(target_protein_ids)
    summaries = []
    mapping_rows = []

    for chunk in pd.read_csv(merge_df_csv, usecols=usecols, dtype=dtype, chunksize=500_000):
        chunk = chunk.loc[chunk["Protein_ID"].isin(target_protein_ids)].copy()
        if chunk.empty:
            continue

        mapping_rows.append(chunk[["virus_id", "Protein_ID"]].drop_duplicates())
        chunk["fimo_scan_length"] = parse_scan_length(chunk["fragment"])
        chunk["context"] = np.where(chunk["chr"].eq("MT"), "Host mitochondrial genome", "Host nuclear genome")
        summaries.append(
            chunk.groupby(["Protein_ID", "context"], as_index=False).agg(
                raw_hit_count=("hit_count", "sum"),
                scan_length_bp=("fimo_scan_length", "sum"),
            )
        )

    if not mapping_rows:
        raise ValueError("No Similar_list Protein_ID values were found in merge_df.csv.")

    mapping = pd.concat(mapping_rows, ignore_index=True).drop_duplicates()
    duplicated = mapping.groupby("virus_id")["Protein_ID"].nunique()
    if (duplicated > 1).any():
        raise ValueError(f"virus_id maps to multiple Protein_ID values: {duplicated[duplicated > 1].to_dict()}")

    host_summary = pd.concat(summaries, ignore_index=True)
    host_summary = (
        host_summary.groupby(["Protein_ID", "context"], as_index=False)
        .agg(raw_hit_count=("raw_hit_count", "sum"), scan_length_bp=("scan_length_bp", "sum"))
    )
    return mapping, host_summary


def fimo_data_row_count(path):
    count = 0
    with path.open() as handle:
        next(handle, None)
        for line in handle:
            if line and not line.startswith("#"):
                count += 1
    return count


def summarize_fimo_files(fimo_dir, protein_to_virus):
    rows = []
    for path in sorted(fimo_dir.glob("*_fimo.tsv")):
        protein_id = path.name[: -len("_fimo.tsv")]
        virus_id = protein_to_virus.get(protein_id)
        if virus_id is None:
            continue
        rows.append(
            {
                "Protein_ID": protein_id,
                "virus_id": virus_id,
                "context": "Viral genome",
                "raw_hit_count": fimo_data_row_count(path),
                "fimo_file": path.name,
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError("No matched viral FIMO TSV files were found.")
    return (
        out.groupby(["Protein_ID", "virus_id", "context"], as_index=False)
        .agg(raw_hit_count=("raw_hit_count", "sum"), fimo_file=("fimo_file", ";".join))
    )


def build_context_table(binding_df, host_summary, viral_summary):
    viral_lengths = binding_df[["Protein_ID", "Genome_length (bp)"]].rename(
        columns={"Genome_length (bp)": "scan_length_bp"}
    )
    viral = viral_summary.merge(viral_lengths, on="Protein_ID", how="left")

    host = host_summary.copy()
    host["virus_id"] = np.nan
    host["fimo_file"] = np.nan

    context_df = pd.concat(
        [
            viral[["Protein_ID", "virus_id", "context", "raw_hit_count", "scan_length_bp", "fimo_file"]],
            host[["Protein_ID", "virus_id", "context", "raw_hit_count", "scan_length_bp", "fimo_file"]],
        ],
        ignore_index=True,
    )
    context_df["context"] = pd.Categorical(context_df["context"], CONTEXTS, ordered=True)
    context_df["raw_hit_count"] = context_df["raw_hit_count"].astype(int)
    context_df["scan_length_bp"] = context_df["scan_length_bp"].astype(int)
    return context_df.sort_values(["Protein_ID", "context"]).reset_index(drop=True)


def hypergeom_one_vs_rest(k_context, n_context_bp, total_hits, total_bp):
    expected = total_hits * n_context_bp / total_bp
    fold = k_context / expected if expected > 0 else np.nan
    p_over = hypergeom.sf(k_context - 1, int(total_bp), int(total_hits), int(n_context_bp))
    return expected, fold, p_over


def adjust_pvalues(out):
    reject_bh, p_bh, _, _ = multipletests(out["hypergeom_overrepresentation_p_raw"], alpha=ALPHA, method="fdr_bh")
    reject_holm, p_holm, _, _ = multipletests(out["hypergeom_overrepresentation_p_raw"], alpha=ALPHA, method="holm")
    out["hypergeom_overrepresentation_p_BH_FDR"] = p_bh
    out["significant_BH_FDR_0.05"] = reject_bh
    out["hypergeom_overrepresentation_p_Holm"] = p_holm
    out["significant_Holm_0.05"] = reject_holm
    return out


def per_protein_tests(context_df):
    rows = []
    for protein_id, sub in context_df.groupby("Protein_ID", observed=True):
        total_hits = int(sub["raw_hit_count"].sum())
        total_bp = int(sub["scan_length_bp"].sum())
        for _, row in sub.iterrows():
            expected, fold, p_over = hypergeom_one_vs_rest(
                int(row["raw_hit_count"]),
                int(row["scan_length_bp"]),
                total_hits,
                total_bp,
            )
            rows.append(
                {
                    "Protein_ID": protein_id,
                    "Genome_length_group": "Similar_list",
                    "context": row["context"],
                    "context_raw_hit_count": int(row["raw_hit_count"]),
                    "context_scan_length_bp": int(row["scan_length_bp"]),
                    "total_raw_hit_count_three_contexts": total_hits,
                    "total_scan_length_bp_three_contexts": total_bp,
                    "expected_raw_hit_count": expected,
                    "fold_enrichment": fold,
                    "hypergeom_overrepresentation_p_raw": p_over,
                }
            )
    return adjust_pvalues(pd.DataFrame(rows))


def similar_list_pooled_tests(context_df):
    pooled = (
        context_df.groupby("context", observed=True, as_index=False)
        .agg(
            context_raw_hit_count=("raw_hit_count", "sum"),
            context_scan_length_bp=("scan_length_bp", "sum"),
        )
        .sort_values("context")
    )
    found_n = int(context_df["Protein_ID"].nunique())
    total_hits = int(pooled["context_raw_hit_count"].sum())
    total_bp = int(pooled["context_scan_length_bp"].sum())

    rows = []
    for _, row in pooled.iterrows():
        expected, fold, p_over = hypergeom_one_vs_rest(
            int(row["context_raw_hit_count"]),
            int(row["context_scan_length_bp"]),
            total_hits,
            total_bp,
        )
        rows.append(
            {
                "Similar_list_found_n": found_n,
                "context": row["context"],
                "context_raw_hit_count": int(row["context_raw_hit_count"]),
                "context_scan_length_bp": int(row["context_scan_length_bp"]),
                "total_raw_hit_count_three_contexts": total_hits,
                "total_scan_length_bp_three_contexts": total_bp,
                "expected_raw_hit_count": expected,
                "fold_enrichment": fold,
                "hypergeom_overrepresentation_p_raw": p_over,
            }
        )
    return adjust_pvalues(pd.DataFrame(rows))


def significance_label(p_value):
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "ns"


def count_label(value):
    if value >= 1e9:
        return f"{value / 1e9:.2f}B"
    if value >= 1e6:
        return f"{value / 1e6:.2f}M"
    if value >= 1e3:
        return f"{value / 1e3:.1f}K"
    return f"{value:.0f}"


def save_pooled_enrichment_plot(pooled, output_path, missing_ids):
    plot_df = pooled.copy()
    plot_df["context"] = pd.Categorical(plot_df["context"], CONTEXTS, ordered=True)
    plot_df = plot_df.sort_values("context")

    sns.set_theme(style="whitegrid", context="notebook")
    fig, ax = plt.subplots(figsize=(7.2, 5.15))
    x = np.arange(len(CONTEXTS))
    y_max = max(1.75, plot_df["fold_enrichment"].max() * 1.28)

    bars = ax.bar(
        x,
        plot_df["fold_enrichment"],
        width=0.58,
        color=[COLORS[context] for context in plot_df["context"].astype(str)],
        edgecolor="black",
        linewidth=0.8,
    )
    for bar, (_, row) in zip(bars, plot_df.iterrows()):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + y_max * 0.025,
            f"{row['fold_enrichment']:.2f}\n{significance_label(row['hypergeom_overrepresentation_p_BH_FDR'])}",
            ha="center",
            va="bottom",
            fontsize=11,
            linespacing=0.95,
        )

    ax.axhline(1, color="black", linestyle="--", linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels([CONTEXT_LABELS[context] for context in CONTEXTS])
    ax.set_ylim(0, y_max)
    ax.set_xlabel("Genomic context")
    ax.set_ylabel("Hypergeometric fold enrichment\n(observed / expected raw FIMO hits)")
    ax.set_title(f"Similar_list pooled hypergeometric enrichment (n={int(plot_df['Similar_list_found_n'].iloc[0])} proteins)")
    ax.text(
        0.5,
        -0.23,
        "*, **, *** indicate BH-FDR adjusted p < 0.05, 0.01, 0.001.\n"
        f"Missing from Binding/raw-count results: {', '.join(missing_ids) if missing_ids else 'none'}",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    fig.savefig(output_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_pooled_observed_expected_heatmap(pooled, output_path, missing_ids):
    plot_df = pooled.copy()
    plot_df["context"] = pd.Categorical(plot_df["context"], CONTEXTS, ordered=True)
    plot_df = plot_df.sort_values("context")
    labels = [CONTEXT_LABELS[context] for context in plot_df["context"].astype(str)]

    observed = pd.DataFrame({"Similar_list pooled": np.log10(plot_df["context_raw_hit_count"].to_numpy() + 1)}, index=labels)
    expected = pd.DataFrame({"Similar_list pooled": np.log10(plot_df["expected_raw_hit_count"].to_numpy() + 1)}, index=labels)
    fold = pd.DataFrame({"Similar_list pooled": plot_df["fold_enrichment"].to_numpy()}, index=labels)
    observed_labels = np.array([[count_label(v)] for v in plot_df["context_raw_hit_count"]])
    expected_labels = np.array([[count_label(v)] for v in plot_df["expected_raw_hit_count"]])
    fold_labels = np.array(
        [[f"{row.fold_enrichment:.2f}{significance_label(row.hypergeom_overrepresentation_p_BH_FDR)}"] for row in plot_df.itertuples()]
    )

    sns.set_theme(style="whitegrid", context="notebook")
    fig, axes = plt.subplots(1, 3, figsize=(10.4, 4.6), sharey=True)

    sns.heatmap(
        observed,
        annot=observed_labels,
        fmt="",
        cmap="Blues",
        linewidths=0.8,
        linecolor="white",
        cbar_kws={"label": "log10(observed raw hits + 1)"},
        ax=axes[0],
    )
    axes[0].set_title("Observed raw hits")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Genomic context")

    sns.heatmap(
        expected,
        annot=expected_labels,
        fmt="",
        cmap="Greens",
        linewidths=0.8,
        linecolor="white",
        cbar_kws={"label": "log10(expected hits + 1)"},
        ax=axes[1],
    )
    axes[1].set_title("Expected hits")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("")

    sns.heatmap(
        fold,
        annot=fold_labels,
        fmt="",
        cmap="YlOrRd",
        vmin=0,
        vmax=max(2.0, float(fold.to_numpy().max())),
        linewidths=0.8,
        linecolor="white",
        cbar_kws={"label": "Fold enrichment"},
        ax=axes[2],
    )
    axes[2].set_title("Fold enrichment")
    axes[2].set_xlabel("")
    axes[2].set_ylabel("")

    fig.suptitle(f"Similar_list pooled observed, expected, and enriched FIMO hits (n={int(plot_df['Similar_list_found_n'].iloc[0])})", y=1.02)
    fig.text(
        0.5,
        -0.02,
        f"Missing from Binding/raw-count results: {', '.join(missing_ids) if missing_ids else 'none'}",
        ha="center",
        va="top",
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    fig.savefig(output_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    binding_csv = args.binding_csv.resolve()
    merge_df_csv = args.merge_df.resolve()
    fimo_dir = args.fimo_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    validate_inputs(binding_csv, merge_df_csv, fimo_dir)

    binding_df = pd.read_csv(binding_csv)
    binding_ids = set(binding_df["Protein_ID"])
    requested_in_binding = [pid for pid in SIMILAR_LIST if pid in binding_ids]
    missing_from_binding = [pid for pid in SIMILAR_LIST if pid not in binding_ids]

    mapping, host_summary = summarize_merge_df(merge_df_csv, requested_in_binding)
    found_ids = [pid for pid in requested_in_binding if pid in set(mapping["Protein_ID"])]
    missing_ids = [pid for pid in SIMILAR_LIST if pid not in set(found_ids)]

    similar_binding_df = binding_df.loc[binding_df["Protein_ID"].isin(found_ids)].copy()
    similar_binding_df["Protein_ID"] = pd.Categorical(similar_binding_df["Protein_ID"], categories=SIMILAR_LIST, ordered=True)
    similar_binding_df = similar_binding_df.sort_values("Protein_ID")
    similar_host_summary = host_summary.loc[host_summary["Protein_ID"].isin(found_ids)].copy()
    viral_summary = summarize_fimo_files(fimo_dir, dict(zip(mapping["Protein_ID"], mapping["virus_id"])))
    viral_summary = viral_summary.loc[viral_summary["Protein_ID"].isin(found_ids)].copy()
    context_df = build_context_table(similar_binding_df, similar_host_summary, viral_summary)

    counts = context_df.groupby("Protein_ID", observed=True)["context"].nunique()
    incomplete_ids = counts[counts != len(CONTEXTS)].index.astype(str).tolist()
    if incomplete_ids:
        raise ValueError(f"Expected three contexts per Similar_list protein; incomplete IDs: {incomplete_ids}")

    per_protein = per_protein_tests(context_df)
    pooled = similar_list_pooled_tests(context_df)

    per_protein_path = output_dir / "Similar_list_three_context_hypergeometric_per_protein_results.csv"
    pooled_path = output_dir / "Similar_list_three_context_hypergeometric_pooled_results.csv"
    found_path = output_dir / "Similar_list_three_context_pooled_found_ids.csv"
    missing_path = output_dir / "Similar_list_three_context_missing_ids.csv"
    enrichment_plot_path = output_dir / "Similar_list_three_context_hypergeometric_pooled_enrichment.pdf"
    heatmap_path = output_dir / "Similar_list_three_context_hypergeometric_pooled_observed_expected_heatmap.pdf"

    per_protein.to_csv(per_protein_path, index=False)
    pooled.to_csv(pooled_path, index=False)
    pd.DataFrame({"found_Protein_ID": found_ids}).to_csv(found_path, index=False)
    pd.DataFrame({"missing_Protein_ID": missing_ids}).to_csv(missing_path, index=False)
    save_pooled_enrichment_plot(pooled, enrichment_plot_path, missing_ids)
    save_pooled_observed_expected_heatmap(pooled, heatmap_path, missing_ids)

    if missing_from_binding:
        print(f"Missing from Binding.csv: {', '.join(missing_from_binding)}")
    print(f"Found Similar_list proteins with raw-count results: {len(found_ids)}")
    print(f"Saved: {per_protein_path}")
    print(f"Saved: {pooled_path}")
    print(f"Saved: {found_path}")
    print(f"Saved: {missing_path}")
    print(f"Saved: {enrichment_plot_path}")
    print(f"Saved: {heatmap_path}")
    print(pooled.to_string(index=False))


if __name__ == "__main__":
    main()
