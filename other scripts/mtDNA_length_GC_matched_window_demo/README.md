# mtDNA-length and GC-matched window demo

Purpose: demonstrate how nuclear control windows were selected with the same length as human mtDNA and similar GC content.

This folder is a small reproducible demo. It does not run FIMO and does not include protein or motif identifiers.

## Reference Data

Use Ensembl release 110 / GRCh38, not `current_fasta`, so the input FASTA files are fixed.

```bash
mkdir -p input
cd input

# Nuclear genome used for nuclear-window sampling, soft-masked GRCh38 primary assembly
wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz

# mtDNA, soft-masked GRCh38 chromosome MT
wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.MT.fa.gz

# Optional checksum file from the same Ensembl release directory
wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/CHECKSUMS

gunzip Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna_sm.chromosome.MT.fa.gz

# Required FASTA index for random access
samtools faidx Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
cd ..
```

If using the unmasked mtDNA FASTA instead, the release-110 file is:

```bash
wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz
```

## Version And Standards

- Nuclear genome version: Homo sapiens GRCh38 primary assembly, Ensembl release 110, soft-masked sequence (`dna_sm`)
- Nuclear sampling scope: chromosomes `1-22`, `X`, and `Y`; chromosome `MT` is excluded from the nuclear control set
- mtDNA version: GRCh38 mitochondrial chromosome `MT`, Ensembl release 110, positions `1-16569`
- Window length standard: all nuclear control windows are `16,569 bp`, equal to the GRCh38 mtDNA length
- Coordinate standard in CSV output: `start_0based` and `end_0based` use 0-based half-open coordinates, `[start, end)`
- Matching standard: selected nuclear windows are ranked by absolute GC difference from mtDNA after excluding windows with excessive `N` content

## Local Paths Used In This Workspace

The analysis in this workspace used these local files:

- Nuclear genome FASTA: `/home/tangmeifang/Human_SELEX/genome/human/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa`
- mtDNA FASTA: `/home/tangmeifang/Human_SELEX/binding/split_human_fasta/NC001538VP2_protein_YP_717937.1_window_16569/MT_dna_sm_chromosome_chromosome_GRCh38_MT_1_16569_1_REF_1_16569.txt`

These local paths are not required for reuse; use `--genome-fasta` and `--mtdna-fasta` to provide your own downloaded files.

## Run

```bash
python run_mtDNA_length_GC_matched_window_demo.py \
  --genome-fasta input/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
  --mtdna-fasta input/Homo_sapiens.GRCh38.dna_sm.chromosome.MT.fa \
  --n-windows 50 \
  --n-candidates 5000 \
  --seed 20260712
```

Outputs are written to `results/`:

- `mtDNA_length_GC_matched_windows_demo.csv`
- `mtDNA_length_GC_matched_windows_demo.fa`
- `mtDNA_length_GC_matched_windows_demo_summary.txt`

Columns in the CSV:

- `window_name`: selected control window name
- `rank_by_gc_match`: rank after sorting candidate windows by GC difference
- `chrom`, `start_0based`, `end_0based`: nuclear genome coordinates
- `length_bp`: window length, equal to mtDNA length
- `gc_fraction`: GC content of the selected nuclear window
- `mtDNA_gc_fraction`: GC content of the mtDNA sequence
- `absolute_gc_difference`: absolute GC difference between the nuclear window and mtDNA
- `n_fraction`: fraction of `N` bases in the selected nuclear window
