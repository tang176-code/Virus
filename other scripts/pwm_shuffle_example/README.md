# Minimal PWM base-label shuffle example

This folder is the concise prepublication version intended to show reviewers how the matched control PWMs were generated.

It contains only a toy example of the base-label-shuffling procedure. It does not include the subsequent FIMO scanning, full mtDNA hit-count data, or statistical analysis.

## File

- `example_pwm_base_label_shuffle.py`: minimal demonstration of the base-label-shuffling procedure used to generate matched control PWMs.

## Run

```bash
python example_pwm_base_label_shuffle.py
```

## Method

For each PWM position, the four probability values are kept unchanged, but the A/C/G/T labels are randomly permuted. This preserves motif length and per-position information content, while disrupting the original nucleotide preference.
