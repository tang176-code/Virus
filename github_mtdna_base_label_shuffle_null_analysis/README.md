# mtDNA base-label PWM shuffle null analysis

This folder contains the source data, compact analysis code, plotting code, and final figure for the mtDNA base-label PWM shuffle null analysis.

The GitHub-facing CSV files omit the internal `ID` column and use `Protein_ID` as the public protein identifier.

## Files

- `example_pwm_base_label_shuffle.py`: minimal toy example showing how A/C/G/T labels were shuffled within each PWM column.
- `mtdna_base_label_shuffle_null_observed_counts.csv`: observed mtDNA FIMO hit counts from the original viral PWMs, after excluding three non-human-infecting virus proteins.
- `mtdna_base_label_shuffle_null_shuffled_counts.csv`: mtDNA hit counts from 1,000 matched base-label-shuffled PWMs per protein.
- `analyze_mtdna_base_label_shuffle_null.py`: compact analysis script that calculates null mean, 95% null interval, fold enrichment, z-score, empirical P value, and BH-FDR.
- `mtdna_base_label_shuffle_null_source_data.csv`: source data used for plotting, with column descriptions at the top of the CSV.
- `plot_mtdna_base_label_shuffle_null_combined.py`: plotting script used to generate the final figure.
- `mtdna_base_label_shuffle_null_combined.pdf`: final publication-style PDF figure.
- `mtdna_base_label_shuffle_null_combined.svg`: editable vector version of the final figure.

## Reproduce

```bash
python analyze_mtdna_base_label_shuffle_null.py
python plot_mtdna_base_label_shuffle_null_combined.py
```

The PWM shuffle null randomly permutes A/C/G/T labels within each PWM column. This preserves motif width and per-position information content, while disrupting the original base preference.
