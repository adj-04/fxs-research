# fxs-research

**Gene regulatory network construction in Fragile X syndrome (FXS)**: an
integrative transcriptomic and epigenomic framework with permutation-tested
networks, validated across two independent patient-derived cohorts.

The pipeline combines RNA-seq differential expression, cell-composition-adjusted
DNA methylation, and four network-inference methods (GENIE3, PANDA, LIONESS,
ALPACA). It tests every network statistic against an FDR-capable permutation
null and validates candidate regulators against ENCODE ChIP-seq, TRRUST and
CollecTRI. Everything runs from public data with one Colab notebook.

Paper: [`paper/fxs-paper-v5.pdf`](paper/fxs-paper-v5.pdf)

---

## Cohorts

| Role | Accession | Material | Design |
|---|---|---|---|
| Primary | [GSE146339](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146339) | iPSC-derived cerebral organoids, day 28 and day 105 | 6 FXS vs 6 control, `~ timepoint + condition` |
| Validation | [GSE218676](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218676) | Patient NPC / B-lymphocyte / iPSC lines | 5 FXS vs 5 control lines (8 donors), `~ cell_type + condition` |
| Methylation | [GSE41273](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41273) | Peripheral blood, Illumina 450K, all male | 9 FXS vs 53 control, `~ condition + cell proportions` |

Both RNA-seq cohorts use NCBI-reprocessed gene-level counts (GRCh38.p13), so
quantification differences cannot explain disagreement between them. The three
datasets come from three different tissues and are used deliberately as three
independent windows onto the same lesion.

---

## Key results

| | Result |
|---|---|
| **FMR1 silencing** | Replicates in all three datasets: log2FC −10.47 (organoids), −8.91 (validation lines), Δβ = 0.58 across 15 promoter probes (blood, cell-adjusted) |
| **Significant differential module** | ALPACA module 4: TFs **BARHL1, OTX1** with targets **NEFL, MPZL2, TAF4B**. FDR = 0.0099; still significant with half the permutations (FDR = 0.03) |
| **Best-validated regulator** | **REST**: rank 1 in ChEA3's ENCODE ChIP-seq library for both GENIE3 and PANDA targets; curated-target enrichment in TRRUST (FDR = 0.004) and CollecTRI (FDR = 5 × 10⁻⁶); top-30 GENIE3 hub |
| **TF activity** | SP1 activity lower in FXS in both cohorts (validation FDR = 0.003); CTNNB1 is the only TF significant in both |
| **Tissue generality** | Only FMR1 and FMR1-AS1 are significant in more than one tissue context |
| **Local centrality** | SP1, TFAP2A, NRF1, REST and FOXG1 sit in the top decile of their own module without being genome-wide hubs |

**Negative results, reported in full.** At n = 6 vs 6, no individual edge is
distinguishable from chance:

| Test | Result |
|---|---|
| GENIE3 edges (pooled null, 20 permutations) | 0 of 8.26 M at FDR < 0.05 |
| PANDA edges (randomized-prior null) | 0 at FDR; 5.2% above null q95, so scores are driven by the priors |
| LIONESS edges (limma, timepoint-adjusted) | 0 of 2.9 M at FDR < 0.05 |
| LIONESS genome-wide shift | 52.4% of edges stronger in FXS, 95% CI [48.3, 53.4] |
| Cross-cohort effect sizes | Pearson r = 0.031 (shrunken: 0.008); sign concordance 50.5% |
| Cross-cohort TF activity | Spearman ρ = −0.02 across 737 TFs |
| Validation DEGs, donor-blocked | 63 at line level, 1 (FMR1) once donor is modelled |

All headline numbers: [`results/v5_headline_numbers.csv`](results/v5_headline_numbers.csv).

---

## Running it

### Google Colab (recommended)

Open [`colab/FXS_GRN_v5_Colab.ipynb`](colab/FXS_GRN_v5_Colab.ipynb) in Colab
with the default **Python 3** runtime and run the cells top to bottom.

- Installs R 4.6 / Bioconductor 3.23 as binaries (via [r2u](https://github.com/eddelbuettel/r2u)).
- Downloads every input from GEO, STRING and TRRUST, so there is nothing to upload.
- Writes data, results, figures, checkpoints and logs to Google Drive.
- Permutation steps checkpoint after every iteration. After a disconnect, rerun
  the setup cells and the same stage; it resumes where it stopped.
- Needs about 6 GB of Drive space. Total run time on the free tier is roughly 4–5 hours.

### Locally

```bash
export FXS_REPO=$PWD FXS_ROOT=$PWD      # FXS_ROOT = where outputs are written
Rscript scripts/01_fetch_data.R         # download all inputs once (MD5 manifest)
Rscript scripts/part1_steps1-8.R        # DE, methylation, GENIE3, PANDA
Rscript scripts/part2_steps9-16.R       # networks, LIONESS, ALPACA, enrichment, validation cohort
Rscript scripts/part3_steps17-19.R      # GENIE3/ALPACA nulls, centrality, ENCODE/TRRUST
Rscript scripts/17b-run.R               # PANDA null
Rscript scripts/part4_steps20-22.R      # tissue generality, CollecTRI, headline numbers
```

Each part runs in a fresh R process, so memory is released between stages
(peak about 6 GB). Settings live in `scripts/00_config.R` and can be overridden
with environment variables, for example `FXS_N_PERM_ALPACA=50` or `FXS_CORES=4`.

**FMR1 gate.** Step 4 stops the pipeline if FMR1 is not significantly
downregulated in FXS, which checks that the sample labels are correct.

### Cytoscape figures (local)

Open Cytoscape first (cap its memory with `-Xmx3g` in `Cytoscape.vmoptions`), then in R:

```r
setwd("~/fxs-research")
source("scripts/14_cytoscape_local.R")
```

This builds the styled ALPACA network via RCy3: the FDR-significant module in
full (gold border) plus the strongest remaining edges. It exports PNG/PDF,
per-module figures, XGMML, the style file and a `.cys` session to
`figures/ALPACA/` and `networks/`.

---

## Pipeline

| Steps | Script | Contents |
|---|---|---|
| 0 | `01_fetch_data.R` | GEO counts and 450K betas, STRING v12 links and protein info, TRRUST; MD5 manifest |
| 1–4 | part 1 | Sample annotation, VST, DESeq2 with apeglm shrinkage, FMR1 gate |
| 5 | part 1 | EpiDISH cell proportions, limma DMPs (adjusted and unadjusted), promoter beta |
| 6 | part 1 | TF list (1,330), timepoint removed for networks, label-free target universe |
| 7 | part 1 | GENIE3: 1,330 TFs → 6,333 genes (5,003 non-TF targets) |
| 8 | part 1 | JASPAR motif prior, STRING PPI, methylation weighting (`1.5 − β`), PANDA |
| 9–10 | part 2 | Network construction, hubs (GENIE3 weighted out-degree, PANDA out-degree) |
| 11 | part 2 | LIONESS (12 networks), per-edge limma `~ timepoint + condition` |
| 12–14 | part 2 | ALPACA modules, GO/KEGG enrichment, network export |
| 15–16 | part 2 | Validation cohort DE with donor-blocked limma-voom; cross-cohort overlap, correlation, RRHO |
| 17 | part 3, 17b | Pooled permutation nulls (GENIE3, PANDA, ALPACA); stratified LIONESS bootstrap |
| 18 | part 3 | Centrality, target-universe coverage, module-local centrality |
| 19 | part 3 | ChEA3 / ENCODE ChIP-seq self-recovery, TRRUST enrichment |
| 20–22 | part 4 | Tissue-generality classification, CollecTRI enrichment and TF activity, headline table |

### Design choices worth knowing

- **Label-free gene selection.** Network targets are the highest-variance
  expressed genes (≥10 counts in ≥3 samples) plus the DEGs, chosen without
  FXS/control labels, so differential tests on the networks are not circular.
- **Timepoint removed before network inference.** Developmental stage dominates
  PC1 (|r| = 0.96) and would otherwise be read as regulation.
- **Pooled permutation nulls.** Comparing each edge only with its own
  permutations gives a smallest possible p of 1/(B+1), far too coarse for FDR
  over millions of edges. Every observed score is instead compared with the
  pooled null of all edges across all permutations (resolution about 10⁻⁸).
- **Methylation is a prior, not matched multi-omics.** The blood cohort shares
  no subjects with the RNA-seq cohorts. Control-group promoter methylation
  enters the PANDA motif prior as a baseline-accessibility weight.

---

## Repository layout

```
colab/      Colab notebook, R installer, package installer
scripts/    00_config.R, 01_fetch_data.R, part1–part4, 17b-run.R, 14_cytoscape_local.R
results/    summary tables (DE, methylation, GENIE3, PANDA, LIONESS, ALPACA, validation, tissue)
figures/    all figures, by analysis
networks/   ALPACA network: XGMML, GraphML, Cytoscape style and session
logs/       full run logs, sessionInfo per stage, package versions
data/       metadata, TF list, MANIFEST.csv (inputs are downloaded, not stored)
paper/      manuscript (LaTeX and PDF)
```

Very large intermediate files (complete GENIE3, PANDA and LIONESS edge lists,
and R checkpoints) are not stored; the scripts regenerate them.

## Environment

R 4.6, Bioconductor 3.23. Exact versions of every package are in
[`logs/package_versions.csv`](logs/package_versions.csv) and
`logs/sessionInfo_*.txt`. Core packages: DESeq2, apeglm, limma, edgeR, minfi,
EpiDISH, GENIE3, pandaR, netZooR, igraph, clusterProfiler, JASPAR2020,
motifmatchr, decoupleR, OmnipathR, RCy3.

netZooR and pandaR contain `class(x) == "ExpressionSet"` comparisons that fail
under R ≥ 4.0; Step 11 patches these at the namespace level before LIONESS runs.

See [`CHANGES.md`](CHANGES.md) for the pipeline changelog.

## Citation

A. Jha and S. Sebastian, "Gene Regulatory Network Construction in Fragile X
Syndrome: An Integrative Transcriptomic, Epigenomic, and Permutation-Tested
Computational Framework Across Two Independent Cohorts," preprint, 2026.
