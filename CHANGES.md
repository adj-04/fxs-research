# Changelog

## v5 (November 2026)

Runs end to end on Google Colab from public data, and tightens the statistics
of every network step. The previous release is preserved under the tag
`v4-paper`.

### Correctness
- **GENIE3 orientation.** `GENIE3()` returns a regulator × target matrix; the
  reshaped edge list now labels regulators and targets accordingly. A check
  stops the run if any regulator is not a TF; per-target weight sums are logged.
- **FDR-capable permutation nulls.** With B permutations, a per-edge empirical
  p cannot fall below 1/(B+1), which is too coarse for BH correction over
  millions of edges. GENIE3, PANDA and ALPACA scores are now compared with the
  pooled null of all edges/nodes across all permutations (constant-memory
  histogram; smallest p about 10⁻⁸). The per-edge statistic is still written in
  the `*_peredge` columns.
- **ALPACA module test** uses 10⁵ random draws (was 10³, which capped the
  smallest module FDR at about 0.25 across 249 modules).
- **STRING mapping** uses STRING's own `protein.info` table (all 16,201 proteins
  mapped; 6,490 TF–TF interactions). The previous `org.Hs.eg.db` ENSEMBLPROT
  route mapped about 5% of high-confidence edges.
- **Enrichment universe** is applied to the relaxed GO/KEGG runs as well as the
  strict ones.

### Analysis
- **Target universe.** GENIE3 targets are all TFs, all DEGs and the top 5,000
  highest-variance non-TF genes passing an expression floor (≥10 normalized
  counts in ≥3 samples): 6,333 genes, 5,003 non-TF targets. PANDA's motif
  target pool comes from the same label-free ranking.
- **Timepoint** is regressed out of the network-inference matrix
  (`removeBatchEffect`, condition protected); LIONESS edges are tested with
  limma `~ timepoint + condition`; the bootstrap is stratified by timepoint.
- **Donor structure** in the validation cohort is handled with limma-voom and
  `duplicateCorrelation(block = donor)` as a sensitivity analysis.
- **Blood cell composition** (EpiDISH, RPC) is included as covariates in the
  methylation model; the unadjusted model is kept as a sensitivity analysis.
- **Effect-size shrinkage** (apeglm) for cross-cohort correlation.
- **Hubs** in GENIE3 are ranked by weighted out-degree.
- **New:** tissue-generality classification with pairwise hypergeometric tests;
  CollecTRI curated-target enrichment; CollecTRI/ULM TF activity in both cohorts;
  a headline-numbers table (`results/v5_headline_numbers.csv`).
- RRHO is computed with the built-in implementation of the original algorithm.

### Reproducibility
- `scripts/01_fetch_data.R` downloads every input (GEO counts and series matrix,
  STRING, TRRUST), validates downloads, and writes `data/MANIFEST.csv` (sizes
  and MD5). If NCBI's annotation table is unavailable, gene symbols are built
  from `org.Hs.eg.db` and the source is recorded in `data/ANNOTATION_SOURCE.txt`.
- `scripts/00_config.R` centralizes paths and settings; all analysis switches
  are environment-variable overridable.
- Every permutation loop checkpoints after each iteration and resumes after
  interruption.
- `colab/` contains the notebook, an r2u-based R 4.6 / Bioconductor 3.23
  installer and a package installer that logs exact versions.
- Cytoscape runs locally via `scripts/14_cytoscape_local.R` (RCy3); Colab writes
  a fallback XGMML/GraphML/PNG.

### Settings to reproduce the previous analysis choices
```
FXS_TARGET_SELECTION=pvalue FXS_ADJUST_TIMEPOINT=FALSE FXS_METH_ADJUST_CELLTYPE=FALSE
```
