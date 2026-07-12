#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Minimal example of the PWM base-label shuffle used in this analysis.

This script provides a minimal demonstration of the base-label-shuffling
procedure used to generate matched control PWMs. It does not include the
subsequent FIMO scanning or statistical analysis.

For each PWM position, the four probability values are kept unchanged, but the
A/C/G/T labels are randomly permuted. This preserves motif length and the
per-position information content, while disrupting the original nucleotide
preference.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

BASES = ["A", "C", "G", "T"]


def information_content(row: np.ndarray) -> float:
    """Information content in bits for one PWM position."""
    p = row[row > 0]
    entropy = -np.sum(p * np.log2(p))
    return 2.0 - entropy


def base_label_shuffle_pwm(pwm: pd.DataFrame, seed: int | None = None) -> pd.DataFrame:
    """Shuffle A/C/G/T labels within each PWM position.

    Input rows are motif positions and columns are A/C/G/T probabilities.
    The set of four probability values in each row is unchanged; only the base
    labels assigned to those values are randomized.
    """
    rng = np.random.default_rng(seed)
    shuffled_rows = []

    for _, row in pwm[BASES].iterrows():
        values = row.to_numpy(dtype=float)
        shuffled_values = values[rng.permutation(len(BASES))]
        shuffled_rows.append(shuffled_values)

    return pd.DataFrame(shuffled_rows, columns=BASES, index=pwm.index)


def main() -> None:
    # Toy PWM: rows are motif positions, columns are A/C/G/T probabilities.
    original_pwm = pd.DataFrame(
        [
            [0.80, 0.10, 0.05, 0.05],
            [0.25, 0.25, 0.25, 0.25],
            [0.05, 0.75, 0.10, 0.10],
            [0.10, 0.10, 0.70, 0.10],
        ],
        columns=BASES,
        index=["pos1", "pos2", "pos3", "pos4"],
    )

    shuffled_pwm = base_label_shuffle_pwm(original_pwm, seed=1)

    print("Original PWM")
    print(original_pwm)
    print("\nBase-label-shuffled PWM")
    print(shuffled_pwm)

    print("\nPer-position information content is preserved:")
    for pos in original_pwm.index:
        original_ic = information_content(original_pwm.loc[pos].to_numpy(dtype=float))
        shuffled_ic = information_content(shuffled_pwm.loc[pos].to_numpy(dtype=float))
        print(f"{pos}: original={original_ic:.6f}, shuffled={shuffled_ic:.6f}")


if __name__ == "__main__":
    main()
