# v5 changes — what changed, why, and what to update in the paper

v5 makes the pipeline run on Google Colab with every input fetched from source,
and removes every limitation that was a property of the **analysis** rather than of
the **data**. What remains is sample size (6 vs 6, 5 vs 5) and the fact that the three
public cohorts share no subjects or tissue.

Every change is switchable through environment variables (see `scripts/00_config.R`),
so the v4 numbers can be regenerated for a sensitivity table.

---

## 1. Correctness fixes (results will change)

### 1.1 GENIE3 regulator/target columns were swapped (Step 7)
`GENIE3()` returns `weightMatrix[regulator, target]`. v4 named `melt()`'s first column
`Target` and second `TF`, so every GENIE3 "TF" was a target and vice versa. The v4
output shows it: non-TF genes (FMR1, FMR1-AS1, TSHR, RERG, SHH, CRABP1 …) have GENIE3
out-degree 82–106 in `GENIE3_centrality.csv`, and `Sum_Weight = 1` for every "TF"
(GENIE3 normalises per **target**). The README caveat about `GENIE3_weighted_degree.csv`
was a symptom of this.

**Affected in the paper:** GENIE3 hub ranking (Fig. 7 left), GENIE3-derived ChEA3/ENCODE
ranks (e.g. SP1 rank 5), all three TRRUST FDR-significant results (SP1 4.3×, HOXA9,
HOXA10 — all GENIE3-derived), the GENIE3 null-filtered network, and Sec. IV-D/VII
wording about GENIE3. PANDA, LIONESS and ALPACA results are not affected
(REST/NRF1 ENCODE rank 1 used PANDA targets).

### 1.2 Permutation FDRs could not pass by construction (Step 17)
v4 gave each edge `p = (#permutations with null ≥ observed + 1) / (N + 1)`:

| Test | Permutations | Smallest possible p | Tests | Could BH ever pass? |
|---|---|---|---|---|
| GENIE3 | 50 | 0.0196 | ~1.7 M edges | No |
| PANDA | 20 | 0.048 | ~340 k edges | No |
| ALPACA | 50 | 0.0196 | thousands of nodes | No |

So "zero edges/nodes survived FDR < 0.05" was certain before the data were analysed —
it is not evidence of absent signal. v5 compares each observed score with the
**pooled** null (all edges × all permutations; resolution ~1e-8), accumulated in a
constant-memory histogram. The per-edge v4 statistic is still computed and saved next
to it. GENIE3 nulls use 250 trees vs 1000 observed, which makes the null wider and the
test conservative.

### 1.3 Enrichment universe (Step 13)
The relaxed GO/KEGG runs had no `universe`, so terms were tested against the whole
genome. Fixed (same background as the strict runs).

---

## 2. Limitations removed

| Paper limitation (Sec. VII) | v5 change | Script / step |
|---|---|---|
| **GENIE3 target coverage** — 1,330/1,340 genes were TFs; TF–TF map; ~0.3 % of regulons eligible | Targets = all TFs + DEGs + top 5,000 most variable non-TF genes (label-independent). PANDA's motif pool is drawn from the same universe instead of the top-3,000 by DE p-value (which selected on the FXS contrast) | part1 Step 7, 8B |
| **Unmodeled covariates (timepoint)** — LIONESS cannot adjust | Timepoint regressed out of the network-inference matrix (`removeBatchEffect`, condition protected); LIONESS edges then tested per edge with limma `~ timepoint + condition`; bootstrap stratified by timepoint | part1 Step 6; part2 Step 11; part3 17D |
| **Unmodeled covariates (donor)** — two donors contribute two lines | limma-voom with `duplicateCorrelation(block = donor)` on the validation cohort; agreement with DESeq2 reported; tissue classification rerun with it | part2 Step 15; part4 Step 20 |
| **Blood methylation without cell-composition adjustment** (not listed, but a standard reviewer request) | EpiDISH RPC proportions as covariates in limma; unadjusted model kept as sensitivity; FMR1 probe comparison table. Cohort is all-male, so X-inactivation cannot confound FMR1 — now stated in the log | part1 Step 5A |
| **Weak genome-wide correlation** partly a noise artefact | apeglm-shrunken LFCs for both cohorts; correlation/sign concordance reported on raw and shrunken LFCs | part1 Step 4; part2 Step 16 |
| **In silico validation scope** — two sources | Adds CollecTRI curated regulons (third source) and CollecTRI/ULM TF-activity inference in **both** cohorts, giving a TF-level cross-cohort replication test that does not depend on individual genes replicating | part4 Step 21 |
| Tissue-generality classification not in the code | Implemented, with pairwise chance expectations (hypergeometric) | part4 Step 20 |
| Cytoscape cannot run on Colab | Colab writes node/edge tables plus a fallback XGMML/PNG; `scripts/14_cytoscape_local.R` then drives Cytoscape on your laptop via RCy3 (style, layout, per-module figures, XGMML, style .xml, .cys session) | part2 Step 14; scripts/14_cytoscape_local.R |
| Manual downloads, hard-coded `setwd()` | `01_fetch_data.R` pulls GEO counts, the GSE41273 beta matrix (GEOquery), STRING, TRRUST; writes `data/MANIFEST.csv` (sizes + md5) | scripts/01_fetch_data.R |

## 3. What still is a limitation (data, not analysis)
1. **Sample size** — 6 vs 6 organoids, 5 vs 5 lines; every network test is power-limited.
2. **No matched multi-omics** — methylation is from an independent blood cohort, so it
   can only enter as a population-level prior.
3. **Tissue differences between cohorts** — now handled as a designed comparison, but a
   "tissue-specific" call can still be under-detection elsewhere.
4. **In silico ≠ experimental** — CollecTRI/ENCODE/TRRUST support plausibility, not causality.

---

## 4. Paper edits after the v5 run
Numbers come from `results/v5_headline_numbers.csv` and the files named below.

- **Abstract / V-C / VI / IX**: replace "no GENIE3 edge, PANDA edge or ALPACA node
  survived FDR < 0.05" with the pooled-null counts (`GENIE3_null_summary.csv`,
  `PANDA_null_summary.csv`, `ALPACA_results_with_significance.csv`). Whatever they show,
  add one sentence in Methods (IV-I) explaining the pooled null and why per-edge nulls
  with 20–50 permutations cannot reach BH significance.
- **IV-D, VII (GENIE3 target coverage)**: delete the TF–TF caveat; report genes/targets.
- **V-B, Fig. 7, V-H, Fig. 20**: regenerate; drop SP1/HOXA TRRUST claims unless they
  reappear with correct orientation.
- **IV-C, V-A**: state the cell-composition-adjusted model and the all-male cohort;
  update DMP counts; cite `DMP_model_comparison.csv` as sensitivity.
- **IV-G, V-D, VII (covariates)**: timepoint-adjusted LIONESS + per-edge limma result;
  donor-blocked validation as sensitivity.
- **IV-L, V-E**: add shrunken-LFC correlation next to the raw one.
- **New results paragraph**: CollecTRI TF-activity cross-cohort correlation and the
  status of SP1/REST/NRF1 activity (`TF_activity_CollecTRI_both_cohorts.csv`).
- **Data and Code Availability**: mention `01_fetch_data.R`, `MANIFEST.csv`, the Colab
  notebook, and `logs/package_versions.csv` (R 4.6 / Bioconductor 3.23).

## 5. Reproducing v4 numbers
```
FXS_TARGET_SELECTION=pvalue FXS_ADJUST_TIMEPOINT=FALSE FXS_METH_ADJUST_CELLTYPE=FALSE
```
(The GENIE3 orientation fix and pooled nulls stay on; the per-edge v4 p-values are
still written in the `*_peredge` columns.)
