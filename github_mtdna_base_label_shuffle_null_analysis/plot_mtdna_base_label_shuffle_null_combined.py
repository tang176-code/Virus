#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Draw base-label shuffle null figures with native editable Matplotlib text."""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent
SOURCE_CSV = BASE / "mtdna_base_label_shuffle_null_source_data.csv"
ALL_STEM = BASE / "mtdna_base_label_shuffle_null_all_proteins"
GROUP_STEM = BASE / "mtdna_base_label_shuffle_null_group_summary"
COMBINED_STEM = BASE / "mtdna_base_label_shuffle_null_combined"
SOURCE_DATA_CSV = BASE / "mtdna_base_label_shuffle_null_source_data.csv"

FONT_DIR = Path("/home/tangmeifang/Human_SELEX/text")
FONT_FILES = [
    FONT_DIR / "MyriadPro-Regular.ttf",
    FONT_DIR / "Myriad Pro Bold.ttf",
    FONT_DIR / "MyriadPro-It.ttf",
]

THRESHOLDS = [0.05, 0.01, 0.005, 0.001]
GROUPS = ["Large", "Mid", "Small"]
GROUP_X = {group: i for i, group in enumerate(GROUPS)}
COLORS = {"Large": "#1F77B4", "Mid": "#E69F00", "Small": "#009E73"}
LEGEND_LABELS = {"Large": "Large", "Mid": "Medium", "Small": "Small"}
MODULE_COLOR = "#D00000"
MODULE_1_PROTEINS = {
    "NP_043126.1",
    "YP_001111256.1",
    "YP_008431134.1",
    "YP_009111420.1",
    "YP_009553535.1",
    "YP_717937.1",
    "YP_717939.1",
}


def configure_matplotlib() -> None:
    for font_file in FONT_FILES:
        if not font_file.exists():
            raise FileNotFoundError(f"Required font file not found: {font_file}")
        font_manager.fontManager.addfont(str(font_file))

    matplotlib.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "font.family": "sans-serif",
        "font.sans-serif": ["Myriad Pro"],
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Myriad Pro",
        "mathtext.it": "Myriad Pro:italic",
        "mathtext.bf": "Myriad Pro:bold",
        "mathtext.cal": "Myriad Pro",
        "mathtext.tt": "Myriad Pro",
        "axes.unicode_minus": False,
    })


def threshold_title(threshold: float) -> str:
    return rf"$\it{{P}}$<{threshold:g}"


CSV_COLUMN_DESCRIPTIONS = [
    ("Protein_ID", "Database accession number of the viral nucleic acid-binding protein."),
    ("genome_length", "Length of the corresponding viral genome, in base pairs."),
    ("Genome_length_group", "Viral genome-size category: Large, Medium, or Small."),
    ("PWM_width", "Width of the PWM, corresponding to the number of nucleotide positions in the motif."),
    ("FIMO_p_threshold", "FIMO P-value filtering threshold used to retain motif hits."),
    ("observed_mtDNA_hits", "Number of mtDNA motif hits detected using the original PWM at the indicated filtering threshold."),
    ("shuffled_null_mean", "Mean mtDNA motif-hit count obtained from the matched base-label-shuffled PWMs."),
    ("shuffled_null_SD", "Standard deviation of mtDNA motif-hit counts in the shuffled-PWM null distribution."),
    ("empirical_null_2.5_percentile", "2.5th percentile of the shuffled-PWM null distribution of mtDNA motif-hit counts."),
    ("empirical_null_97.5_percentile", "97.5th percentile of the shuffled-PWM null distribution of mtDNA motif-hit counts."),
    ("fold_enrichment", "Observed mtDNA hit count divided by the shuffled-PWM null mean; values above 1 indicate enrichment relative to the null mean."),
    ("z_score", "Standardized difference, calculated as (observed hits - null mean) / null SD."),
    ("n_shuffled_ge_observed", "Number of shuffled PWMs with mtDNA hit counts greater than or equal to the observed count."),
    ("empirical_P_over", "One-sided empirical P value, calculated as (1 + number of shuffled PWMs with hits >= observed hits) / (1 + total shuffled PWMs)."),
    ("n_shuffled_PWMs", "Number of matched shuffled PWMs generated for each original PWM."),
    ("n_identity_shuffles_rejected", "Number of identity shuffle attempts that were rejected and regenerated."),
    ("BH_adjusted_empirical_P", "BH-adjusted empirical P value, calculated separately across proteins at each FIMO threshold."),
    ("significant_BH_FDR_0.05", "Indicates whether the result remained significant after BH correction at FDR < 0.05."),
]

SOURCE_COLUMN_RENAME = {
    "width": "PWM_width",
    "p_threshold": "FIMO_p_threshold",
    "observed": "observed_mtDNA_hits",
    "null_mean": "shuffled_null_mean",
    "null_sd": "shuffled_null_SD",
    "null_ci_lo": "empirical_null_2.5_percentile",
    "null_ci_hi": "empirical_null_97.5_percentile",
    "n_null_ge_obs": "n_shuffled_ge_observed",
    "emp_p_over": "empirical_P_over",
    "n_shuffle": "n_shuffled_PWMs",
    "n_identity_skipped": "n_identity_shuffles_rejected",
    "emp_p_BH_FDR": "BH_adjusted_empirical_P",
    "significant_BH_0.05": "significant_BH_FDR_0.05",
}

SOURCE_OUTPUT_COLUMNS = [name for name, _description in CSV_COLUMN_DESCRIPTIONS]


def write_source_data_csv(data: pd.DataFrame) -> None:
    source = data.rename(columns=SOURCE_COLUMN_RENAME).copy()
    source["Genome_length_group"] = source["Genome_length_group"].astype(str).replace({"Mid": "Medium"})
    source = source[SOURCE_OUTPUT_COLUMNS]
    with SOURCE_DATA_CSV.open("w", encoding="utf-8", newline="") as handle:
        handle.write("# CSV column descriptions\n")
        handle.write("# Lines beginning with # are metadata comments. To read the data table with pandas, use pd.read_csv(path, comment='#').\n")
        for name, description in CSV_COLUMN_DESCRIPTIONS:
            handle.write(f"# {name}: {description}\n")
        handle.write("#\n")
        source.to_csv(handle, index=False)


def load_data() -> pd.DataFrame:
    data = pd.read_csv(SOURCE_CSV, comment="#").rename(columns={
        "PWM_width": "width",
        "FIMO_p_threshold": "p_threshold",
        "observed_mtDNA_hits": "observed",
        "shuffled_null_mean": "null_mean",
        "shuffled_null_SD": "null_sd",
        "empirical_null_2.5_percentile": "null_ci_lo",
        "empirical_null_97.5_percentile": "null_ci_hi",
        "n_shuffled_ge_observed": "n_null_ge_obs",
        "empirical_P_over": "emp_p_over",
        "n_shuffled_PWMs": "n_shuffle",
        "n_identity_shuffles_rejected": "n_identity_skipped",
        "BH_adjusted_empirical_P": "emp_p_BH_FDR",
        "significant_BH_FDR_0.05": "significant_BH_0.05",
    }).copy()
    data["Genome_length_group"] = data["Genome_length_group"].astype(str).replace({"Medium": "Mid"})
    data["Genome_length_group"] = pd.Categorical(data["Genome_length_group"], GROUPS, ordered=True)
    data["ci_lo_fold"] = data["null_ci_lo"] / data["null_mean"].replace(0, np.nan)
    data["ci_hi_fold"] = data["null_ci_hi"] / data["null_mean"].replace(0, np.nan)
    data["Highlighted_module1"] = np.where(data["Protein_ID"].isin(MODULE_1_PROTEINS), "Yes", "No")
    data = data.sort_values(["p_threshold", "Genome_length_group", "genome_length"], ascending=[False, True, False])
    write_source_data_csv(data)
    return data


def ordered_proteins(data: pd.DataFrame) -> pd.DataFrame:
    first = data[data["p_threshold"].eq(THRESHOLDS[0])].copy()
    first = first.sort_values(["Genome_length_group", "genome_length"], ascending=[True, False])
    return first[["Protein_ID", "Genome_length_group", "genome_length"]].reset_index(drop=True)


def x_max_for_forest(data: pd.DataFrame) -> float:
    x_max = max(2.1, float(data["ci_hi_fold"].max()) * 1.08, float(data["fold_enrichment"].max()) * 1.12)
    return float(np.ceil(x_max * 2) / 2)


def legend_handles(include_median: bool = False):
    handles = [
        mlines.Line2D([], [], color=COLORS[g], marker="o", linestyle="None", markersize=4,
                      markeredgecolor="black", markeredgewidth=0.35, label=LEGEND_LABELS[g])
        for g in GROUPS
    ]
    handles.extend([
        mlines.Line2D([], [], color=MODULE_COLOR, marker="D", linestyle="None", markersize=4,
                      markeredgecolor="black", markeredgewidth=0.35, label="Module 1"),
        mlines.Line2D([], [], color="black", linestyle="--", linewidth=0.5, label="Fold=1"),
    ])
    if include_median:
        handles.append(mlines.Line2D([], [], color="black", linestyle="-", linewidth=0.65, label="Median"))
    return handles


def style_axis(ax) -> None:
    ax.grid(False)
    ax.tick_params(axis="both", which="major", direction="out", length=3.2, width=0.5, colors="black")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(0.5)


def draw_all_proteins_on_axes(data: pd.DataFrame, axes, *, show_title: bool = False) -> None:
    order = ordered_proteins(data)
    protein_order = order["Protein_ID"].tolist()
    y_positions = {protein_id: len(protein_order) - 1 - i for i, protein_id in enumerate(protein_order)}
    y = np.arange(len(protein_order))
    label_rows = order.iloc[::-1].reset_index(drop=True)
    module_ids = set(data.loc[data["Highlighted_module1"].eq("Yes"), "Protein_ID"])
    group_lookup = dict(zip(order["Protein_ID"], order["Genome_length_group"]))
    x_max = x_max_for_forest(data)

    for panel_index, (ax, threshold) in enumerate(zip(axes, THRESHOLDS)):
        panel = data[data["p_threshold"].eq(threshold)].copy()
        panel = panel.set_index("Protein_ID").reindex(protein_order).reset_index()

        for group in GROUPS:
            group_positions = [y_positions[p] for p in order.loc[order["Genome_length_group"].eq(group), "Protein_ID"]]
            ax.axhspan(min(group_positions) - 0.5, max(group_positions) + 0.5,
                       facecolor=COLORS[group], edgecolor=COLORS[group], alpha=0.055,
                       linewidth=0.4, zorder=-10)

        for group in GROUPS:
            sub = panel[panel["Genome_length_group"].eq(group)].copy()
            regular = sub[~sub["Highlighted_module1"].eq("Yes")]
            if len(regular):
                yy = np.array([y_positions[p] for p in regular["Protein_ID"]])
                ax.hlines(yy, regular["ci_lo_fold"], regular["ci_hi_fold"],
                          color=COLORS[group], alpha=0.78, linewidth=0.85, zorder=1)
                ax.scatter(regular["fold_enrichment"], yy, s=16, marker="o", color=COLORS[group],
                           edgecolor="black", linewidth=0.35, zorder=3)

            module = sub[sub["Highlighted_module1"].eq("Yes")]
            if len(module):
                yy = np.array([y_positions[p] for p in module["Protein_ID"]])
                ax.hlines(yy, module["ci_lo_fold"], module["ci_hi_fold"],
                          color=MODULE_COLOR, alpha=0.98, linewidth=1.35, zorder=2)
                ax.scatter(module["fold_enrichment"], yy, s=20, marker="D", color=MODULE_COLOR,
                           edgecolor="black", linewidth=0.35, zorder=4)

        ax.axvline(1.0, color="black", linestyle="--", linewidth=0.5, zorder=0)
        ax.set_title(threshold_title(threshold), pad=5)
        ax.set_xlabel("Fold enrichment")
        ax.set_xlim(0, x_max)
        ax.set_ylim(-0.8, len(protein_order) - 0.2)
        ax.set_xticks(np.arange(0, x_max + 0.001, 1.0))
        ax.set_yticks(y)
        if panel_index == 0:
            ax.set_yticklabels([str(row.Protein_ID) for row in label_rows.itertuples()])
            for tick_label, row in zip(ax.get_yticklabels(), label_rows.itertuples()):
                protein_id = row.Protein_ID
                tick_label.set_color(MODULE_COLOR if protein_id in module_ids else COLORS[group_lookup[protein_id]])
                tick_label.set_fontweight("normal")
            ax.set_ylabel("Protein ID")
        else:
            ax.tick_params(labelleft=False)
            ax.set_ylabel("")
        style_axis(ax)


def jitter_positions(sub: pd.DataFrame) -> np.ndarray:
    x = np.array([GROUP_X[group] for group in sub["Genome_length_group"]], dtype=float)
    offsets = np.zeros(len(sub), dtype=float)
    for group in GROUPS:
        idx = np.where(sub["Genome_length_group"].astype(str).to_numpy() == group)[0]
        n = len(idx)
        if n:
            offsets[idx] = np.linspace(-0.16, 0.16, n) if n > 1 else 0.0
    return x + offsets


def group_y_limits(data: pd.DataFrame) -> tuple[float, float]:
    y_max = max(2.1, float(data["fold_enrichment"].max()) * 1.12)
    y_max = float(np.ceil(y_max * 2) / 2)
    y_min = max(0.0, float(data["fold_enrichment"].min()) * 0.82)
    if y_min < 0.5:
        y_min = 0.2
    return y_min, y_max


def draw_group_summary_on_axes(data: pd.DataFrame, axes) -> None:
    y_min, y_max = group_y_limits(data)

    for panel_index, (ax, threshold) in enumerate(zip(axes, THRESHOLDS)):
        panel = data[data["p_threshold"].eq(threshold)].copy()
        panel = panel.sort_values(["Genome_length_group", "genome_length"], ascending=[True, False]).reset_index(drop=True)

        for group in GROUPS:
            x = GROUP_X[group]
            ax.axvspan(x - 0.5, x + 0.5, facecolor=COLORS[group], edgecolor=COLORS[group],
                       alpha=0.055, linewidth=0.4, zorder=-10)

        for group in GROUPS:
            sub = panel[(panel["Genome_length_group"].eq(group)) & (panel["Highlighted_module1"].eq("No"))]
            if len(sub):
                ax.scatter(jitter_positions(sub), sub["fold_enrichment"], s=16, marker="o",
                           color=COLORS[group], edgecolor="black", linewidth=0.35, alpha=0.96, zorder=3)

        module = panel[panel["Highlighted_module1"].eq("Yes")]
        if len(module):
            ax.scatter(jitter_positions(module), module["fold_enrichment"], s=20, marker="D",
                       color=MODULE_COLOR, edgecolor="black", linewidth=0.35, alpha=0.98, zorder=4)

        median_label_offset = 0.035 * (y_max - y_min)
        for group in GROUPS:
            sub = panel[panel["Genome_length_group"].eq(group)]
            if len(sub):
                median = float(sub["fold_enrichment"].median())
                ax.hlines(median, GROUP_X[group] - 0.24, GROUP_X[group] + 0.24,
                          color="black", linewidth=0.65, zorder=5)
                label_y = median + median_label_offset
                va = "bottom"
                if label_y > y_max - median_label_offset * 0.4:
                    label_y = median - median_label_offset
                    va = "top"
                ax.text(GROUP_X[group], label_y, f"{median:.2f}",
                        ha="center", va=va, fontsize=8, color="black", zorder=6)

        ax.axhline(1.0, color="black", linestyle="--", linewidth=0.5, zorder=1)
        ax.set_title(threshold_title(threshold), pad=5)
        ax.set_xlim(-0.55, 2.55)
        ax.set_ylim(y_min, y_max)
        ax.set_xticks([GROUP_X[g] for g in GROUPS])
        ax.set_xticklabels([LEGEND_LABELS[g] for g in GROUPS])
        ax.set_xlabel("Genome size group")
        if panel_index == 0:
            ax.set_ylabel("Fold enrichment")
        else:
            ax.set_ylabel("")
            ax.tick_params(labelleft=False)
        style_axis(ax)


def draw_all_proteins(data: pd.DataFrame):
    configure_matplotlib()
    n_shuffle = int(data["n_shuffle"].dropna().max())
    fig, axes = plt.subplots(1, 4, figsize=(9.2, 5.25), sharey=True)
    draw_all_proteins_on_axes(data, axes)
    fig.legend(handles=legend_handles(False), title="Genome size", loc="center right",
               bbox_to_anchor=(0.975, 0.52), frameon=False, fontsize=8, title_fontsize=8)
    fig.suptitle(f"PWM shuffle null across FIMO thresholds (N={n_shuffle:,} shuffles)", y=0.975, fontsize=8)
    fig.text(0.44, 0.018,
             "Points show mtDNA observed/null mean; horizontal lines show base-label shuffle null 95% CI scaled by null mean.",
             ha="center", va="bottom", fontsize=8)
    fig.tight_layout(rect=[0.0, 0.06, 0.895, 0.955], w_pad=0.55)
    return fig


def draw_group_summary(data: pd.DataFrame):
    configure_matplotlib()
    fig, axes = plt.subplots(1, 4, figsize=(9.2, 2.25), sharey=True)
    draw_group_summary_on_axes(data, axes)
    fig.legend(handles=legend_handles(True), title="Genome size", loc="center right",
               bbox_to_anchor=(0.985, 0.55), frameon=False, fontsize=8, title_fontsize=8)
    fig.tight_layout(rect=[0.0, 0.0, 0.875, 0.98], w_pad=0.65)
    return fig


def draw_combined(data: pd.DataFrame):
    configure_matplotlib()
    n_shuffle = int(data["n_shuffle"].dropna().max())
    fig = plt.figure(figsize=(8.27, 11.69))

    # Fit the two native panels onto portrait A4 without changing their visual
    # aspect ratio. The drawing block height follows (5.25 + 2.25) / 9.2 of
    # the available block width, leaving page whitespace instead of stretching.
    gs = fig.add_gridspec(
        2, 4,
        left=0.105,
        right=0.855,
        bottom=0.285,
        top=0.715,
        hspace=0.30,
        wspace=0.13,
        height_ratios=[5.25, 2.25],
    )
    top_axes = [fig.add_subplot(gs[0, i]) for i in range(4)]
    bottom_axes = []
    for i in range(4):
        sharey = bottom_axes[0] if bottom_axes else None
        bottom_axes.append(fig.add_subplot(gs[1, i], sharey=sharey))

    draw_all_proteins_on_axes(data, top_axes)
    draw_group_summary_on_axes(data, bottom_axes)

    top_pos = top_axes[0].get_position()
    bottom_pos = bottom_axes[0].get_position()
    top_center_y = (top_pos.y0 + top_pos.y1) / 2
    bottom_center_y = (bottom_pos.y0 + bottom_pos.y1) / 2

    fig.legend(handles=legend_handles(False), title="Genome size", loc="center left",
               bbox_to_anchor=(0.875, top_center_y), frameon=False, fontsize=8, title_fontsize=8)
    fig.legend(handles=legend_handles(True), title="Genome size", loc="center left",
               bbox_to_anchor=(0.875, bottom_center_y), frameon=False, fontsize=8, title_fontsize=8)
    fig.suptitle(f"PWM shuffle null across FIMO thresholds (N={n_shuffle:,} shuffles)",
                 y=top_pos.y1 + 0.028, fontsize=8)
    fig.text(0.43, top_pos.y0 - 0.028,
             "Points show mtDNA observed/null mean; horizontal lines show base-label shuffle null 95% CI scaled by null mean.",
             ha="center", va="bottom", fontsize=8)
    return fig

def save_figure(fig, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".pdf"), dpi=600)
    fig.savefig(stem.with_suffix(".svg"), dpi=600)
    fig.savefig(stem.with_suffix(".png"), dpi=600)
    plt.close(fig)


def main() -> None:
    data = load_data()
    save_figure(draw_all_proteins(data), ALL_STEM)
    save_figure(draw_group_summary(data), GROUP_STEM)
    save_figure(draw_combined(data), COMBINED_STEM)

    print(f"Source data: {SOURCE_DATA_CSV}")
    for stem in (ALL_STEM, GROUP_STEM, COMBINED_STEM):
        print(f"Wrote {stem.with_suffix('.pdf')}")
        print(f"Wrote {stem.with_suffix('.svg')}")
        print(f"Wrote {stem.with_suffix('.png')}")


if __name__ == "__main__":
    main()
