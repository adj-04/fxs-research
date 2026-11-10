# Changelog

## v5.1 (November 2026)

Addresses an external review of the v5 network statistics. New script
`scripts/part5_robust_nulls.R` (resumable; outputs in `results/robust_nulls/`
and `figures/robust_nulls/`).

### Correctness
- **ALPACA orientation.** `pandaToAlpaca()` treats its first network as the
  control (reference). v5 passed FXS first; parts 2 and 3 now pass control first.
- **ALPACA module test replaced.** `alpaca()` normalises node scores to sum to 2
  within each module, so a module's mean score is 2/size and the v5 pooled
  module test made small modules significant by construction (the five-node
  "module 4", FDR 0.0099). Modules are now scored by raw differential modularity
  Q_m and tested by FXS/control label permutation within timepoint (400 possible
  labelings; 100 run), with a max-statistic family-wise p. Result: no
  differential module (global p = 0.81; best module p = 0.38). The v5 test is
  disabled (`FXS_RUN_OLD_ALPACA_NULL=TRUE` reproduces it); its outputs are in
  `results/ALPACA/superseded_v5/`.
- **GENIE3 null calibration.** The v5 null used 250 trees against 1,000 for the
  observed network, which widened it (2% of observed edges at p < 0.05). The
  null now uses 1,000 trees and is checked by leave-one-out (4.97% of held-out
  null edges at p < 0.05). Observed: 5.13% at p < 0.05, 0 at FDR < 0.10.
- **PANDA comparison like for like.** v5 compared positive observed edges with
  all null edges. All observed edges are now compared with all null edges:
  5.17% above the null q95, the fraction expected by chance.

### Analysis
- **New PANDA shuffled-expression null:** real priors, expression permuted per
  gene. Shuffled networks agree with each other (r = 0.998) more than with the
  observed network (r = 0.894), so the expression data do reshape PANDA.
- Module-local centrality is recomputed on the corrected modules as a
  percentile among TFs only (reference TFs: 16th–78th percentile).
- Headline table rows for the v5 ALPACA tests are replaced by the part 5 results.
- The Cytoscape script reads the corrected modules.
- The GO/KEGG enrichment of v5 ALPACA modules (part 2) used the reversed
  orientation and is no longer reported.

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
- **ALPACA module test** uses 10⁵ random draws (was 10³). Superseded in v5.1.
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
