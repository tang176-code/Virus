#!/usr/bin/env python3
"""Demo: sample mtDNA-length nuclear windows matched by GC content.

This script is intentionally limited to the window-selection step. It does not
run FIMO and does not use protein or motif identifiers. The output shows how
control nuclear windows were selected to have:

  1. the same length as human mtDNA
  2. GC content closest to human mtDNA among randomly sampled nuclear windows

Default paths match the local Human_SELEX workspace. Use command-line options to
make the demo portable.
"""

from __future__ import annotations

import argparse
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_GENOME_FASTA = SCRIPT_DIR / "input" / "Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
DEFAULT_MTDNA_FASTA = SCRIPT_DIR / "input" / "Homo_sapiens.GRCh38.dna_sm.chromosome.MT.fa"
DEFAULT_SEED = 20260712


def read_fasta_sequence(path: Path) -> str:
    parts = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            parts.append(line.upper())
    seq = "".join(parts)
    return "".join(base if base in "ACGTN" else "N" for base in seq)


def wrap_sequence(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def gc_fraction(seq: str) -> float:
    acgt = sum(seq.count(base) for base in "ACGT")
    if acgt == 0:
        return float("nan")
    return (seq.count("G") + seq.count("C")) / acgt


def n_fraction(seq: str) -> float:
    return seq.count("N") / len(seq)


def read_fai(fasta: Path) -> OrderedDict[str, dict[str, int]]:
    fai = Path(str(fasta) + ".fai")
    if not fai.exists():
        raise FileNotFoundError(
            f"Missing FASTA index: {fai}\n"
            "Create it first, for example: samtools faidx <genome.fa>"
        )

    records: OrderedDict[str, dict[str, int]] = OrderedDict()
    with fai.open() as handle:
        for line in handle:
            name, length, offset, line_bases, line_width = line.split()[:5]
            records[name] = {
                "length": int(length),
                "offset": int(offset),
                "line_bases": int(line_bases),
                "line_width": int(line_width),
            }
    return records


def fetch_fasta_window(handle, fai: dict[str, dict[str, int]], chrom: str, start0: int, length: int) -> str:
    rec = fai[chrom]
    line_bases = rec["line_bases"]
    line_width = rec["line_width"]
    byte_offset = rec["offset"] + (start0 // line_bases) * line_width + (start0 % line_bases)
    n_bytes = length + (length // line_bases) + 100
    handle.seek(byte_offset)
    chunk = handle.read(n_bytes)
    seq = chunk.replace(b"\n", b"").replace(b"\r", b"")[:length].upper().decode("ascii")
    if len(seq) != length:
        raise RuntimeError(f"Short FASTA fetch for {chrom}:{start0}-{start0 + length}; got {len(seq)} bp")
    return "".join(base if base in "ACGTN" else "N" for base in seq)


def sample_gc_matched_windows(
    genome_fasta: Path,
    mt_seq: str,
    n_windows: int,
    n_candidates: int,
    seed: int,
    max_n_fraction: float,
) -> pd.DataFrame:
    fai = read_fai(genome_fasta)
    window_len = len(mt_seq)
    mt_gc = gc_fraction(mt_seq)

    chromosomes = [str(i) for i in range(1, 23)] + ["X", "Y"]
    chromosomes = [chrom for chrom in chromosomes if chrom in fai and fai[chrom]["length"] >= window_len]
    if not chromosomes:
        raise RuntimeError("No nuclear chromosomes are long enough for mtDNA-length windows.")

    weights = np.array([fai[chrom]["length"] for chrom in chromosomes], dtype=float)
    weights = weights / weights.sum()
    rng = np.random.default_rng(seed)

    candidates = []
    seen = set()
    attempts = 0
    max_attempts = n_candidates * 50

    with genome_fasta.open("rb") as handle:
        while len(candidates) < n_candidates and attempts < max_attempts:
            attempts += 1
            chrom = str(rng.choice(chromosomes, p=weights))
            max_start = fai[chrom]["length"] - window_len
            start0 = int(rng.integers(0, max_start + 1))
            key = (chrom, start0)
            if key in seen:
                continue
            seen.add(key)

            seq = fetch_fasta_window(handle, fai, chrom, start0, window_len)
            nf = n_fraction(seq)
            if nf > max_n_fraction:
                continue
            gc = gc_fraction(seq)
            candidates.append(
                {
                    "chrom": chrom,
                    "start_0based": start0,
                    "end_0based": start0 + window_len,
                    "length_bp": window_len,
                    "gc_fraction": gc,
                    "mtDNA_gc_fraction": mt_gc,
                    "absolute_gc_difference": abs(gc - mt_gc),
                    "n_fraction": nf,
                    "sequence": seq,
                }
            )

    if len(candidates) < n_windows:
        raise RuntimeError(f"Only {len(candidates)} usable candidates were sampled; need {n_windows}.")

    selected = (
        pd.DataFrame(candidates)
        .sort_values(["absolute_gc_difference", "n_fraction", "chrom", "start_0based"])
        .head(n_windows)
        .reset_index(drop=True)
    )
    selected.insert(0, "rank_by_gc_match", np.arange(1, len(selected) + 1))
    selected.insert(1, "window_name", [f"mtDNA_length_GC_match_{i:04d}" for i in range(1, len(selected) + 1)])
    return selected


def write_outputs(selected: pd.DataFrame, mt_seq: str, out_dir: Path, n_candidates: int, seed: int) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / "mtDNA_length_GC_matched_windows_demo.csv"
    fasta_path = out_dir / "mtDNA_length_GC_matched_windows_demo.fa"
    summary_path = out_dir / "mtDNA_length_GC_matched_windows_demo_summary.txt"

    table = selected.drop(columns=["sequence"])
    table.to_csv(csv_path, index=False)

    with fasta_path.open("w") as handle:
        for row in selected.itertuples(index=False):
            header = (
                f">{row.window_name}|{row.chrom}:{row.start_0based}-{row.end_0based}"
                f"|length={row.length_bp}|GC={row.gc_fraction:.6f}"
            )
            handle.write(header + "\n")
            handle.write(wrap_sequence(row.sequence) + "\n")

    with summary_path.open("w") as handle:
        handle.write("mtDNA-length and GC-matched nuclear window demo\n")
        handle.write("=" * 60 + "\n")
        handle.write(f"Random seed: {seed}\n")
        handle.write(f"Candidate nuclear windows sampled: {n_candidates}\n")
        handle.write(f"Selected windows: {len(selected)}\n")
        handle.write(f"Window length: {len(mt_seq)} bp\n")
        handle.write(f"mtDNA GC fraction: {gc_fraction(mt_seq):.6f}\n")
        handle.write(f"Selected nuclear GC median: {selected['gc_fraction'].median():.6f}\n")
        handle.write(
            "Selected nuclear GC range: "
            f"{selected['gc_fraction'].min():.6f}-{selected['gc_fraction'].max():.6f}\n"
        )
        handle.write(
            "Absolute GC difference median: "
            f"{selected['absolute_gc_difference'].median():.6g}\n"
        )
        handle.write("\nOutput files:\n")
        handle.write(f"  {csv_path.name}\n")
        handle.write(f"  {fasta_path.name}\n")

    print(f"Wrote {csv_path}")
    print(f"Wrote {fasta_path}")
    print(f"Wrote {summary_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Sample nuclear windows with mtDNA length and closest GC content."
    )
    parser.add_argument(
        "--genome-fasta",
        type=Path,
        default=DEFAULT_GENOME_FASTA,
        help="GRCh38 nuclear genome FASTA. Default: input/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
    )
    parser.add_argument(
        "--mtdna-fasta",
        type=Path,
        default=DEFAULT_MTDNA_FASTA,
        help="GRCh38 mtDNA FASTA. Default: input/Homo_sapiens.GRCh38.dna_sm.chromosome.MT.fa",
    )
    parser.add_argument("--n-windows", type=int, default=50)
    parser.add_argument("--n-candidates", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--max-n-fraction", type=float, default=0.01)
    parser.add_argument("--out-dir", type=Path, default=Path(__file__).resolve().parent / "results")
    args = parser.parse_args()

    mt_seq = read_fasta_sequence(args.mtdna_fasta)
    selected = sample_gc_matched_windows(
        genome_fasta=args.genome_fasta,
        mt_seq=mt_seq,
        n_windows=args.n_windows,
        n_candidates=args.n_candidates,
        seed=args.seed,
        max_n_fraction=args.max_n_fraction,
    )
    write_outputs(selected, mt_seq, args.out_dir, args.n_candidates, args.seed)


if __name__ == "__main__":
    main()
