# Similar_list three-context pooled hypergeometric enrichment

Self-contained code and data to reproduce the pooled hypergeometric enrichment
analysis of viral protein DNA-binding motif hits (raw FIMO counts) across three
genomic contexts:

1. **Viral genome**
2. **Host mitochondrial genome (mtDNA)**
3. **Host nuclear genome**

For each context, a one-vs-rest hypergeometric over-representation test compares
observed motif hit counts against the expectation from scanned sequence length,
with BH-FDR (primary) and Holm multiple-testing correction. The headline figure
is `results/Similar_list_three_context_hypergeometric_pooled_enrichment.pdf`.

Of the 9 proteins in the Similar_list, 7 have binding/raw-count results and
enter the analysis; `YP_010782993.1` and `AVM80381.1` are absent from
`Binding.csv` and are excluded (reported in
`results/Similar_list_three_context_missing_ids.csv` and in the figure footnote).

## Repository layout

```
├── Similar_list_three_context_hypergeometric_pooled_enrichment.py
├── requirements.txt
├── data/
│   ├── Binding.csv                     # protein list with viral genome lengths
│   ├── merge_df_similar_list.csv.xz    # per-fragment FIMO hit counts on the host
│   │                                   #   genome, subset to the 7 analyzed
│   │                                   #   Similar_list proteins (xz-compressed;
│   │                                   #   pandas reads it directly)
│   └── PWM_fimo_result/
│       └── <Protein_ID>_fimo.tsv       # FIMO scan of each protein's viral genome
└── results/                            # outputs of the script (committed for reference)
    ├── Similar_list_three_context_hypergeometric_pooled_results.csv
    ├── Similar_list_three_context_hypergeometric_per_protein_results.csv
    ├── Similar_list_three_context_pooled_found_ids.csv
    ├── Similar_list_three_context_missing_ids.csv
    ├── Similar_list_three_context_hypergeometric_pooled_enrichment.pdf/.png
    └── Similar_list_three_context_hypergeometric_pooled_observed_expected_heatmap.pdf/.png
```

## Reproduce

```bash
pip install -r requirements.txt   # Python >= 3.9
python Similar_list_three_context_hypergeometric_pooled_enrichment.py
```

This regenerates every file in `results/` from the raw counts in `data/`.
Custom paths can be supplied with `--binding-csv`, `--merge-df`, `--fimo-dir`,
and `--output-dir`.

## Method summary

Pooling the 7 proteins, with total motif hits `K` over total scanned length
`N` bp, the hit count `k` in a context of scanned length `n` bp is tested for
over-representation with `P = scipy.stats.hypergeom.sf(k - 1, N, K, n)`;
expected hits are `K * n / N` and fold enrichment is `k / expected`. P-values
are adjusted with Benjamini-Hochberg FDR (alpha = 0.05, shown as asterisks in
the figure) and Holm methods. The same test is also run per protein
(`results/..._per_protein_results.csv`).

## Data provenance

- `data/merge_df_similar_list.csv.xz` is the subset (columns `fragment`,
  `hit_count`, `chr`, `virus_id`, `Protein_ID`; rows for the 7 analyzed
  proteins) of the full per-fragment host-genome scan table. Fragment IDs
  encode scanned coordinates (`..._REF_<start>_<stop>`), from which scanned
  length is derived. `chr == MT` is assigned to the mtDNA context; all other
  chromosomes to the nuclear context.
- `data/PWM_fimo_result/<Protein_ID>_fimo.tsv` are FIMO scans of each viral
  genome with the corresponding protein's HT-SELEX-derived PWM; viral genome
  lengths for the expected-count calculation come from `data/Binding.csv`.
