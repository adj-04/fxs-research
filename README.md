# fxs-research

Multi-omics regulatory network analysis of Fragile X Syndrome, integrating
RNA-seq differential expression, DNA methylation, and gene regulatory network
inference across two independent patient-derived cohorts.

## Cohorts

| Role | Accession | Material | Design |
|---|---|---|---|
| Primary | [GSE146339](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146339) | iPSC-derived cerebral organoids, day 28 and day 105 | 6 FXS vs 6 control, `~ timepoint + condition` |
| Validation | [GSE218676](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218676) | Patient NPC / B-lymphocyte / iPSC | 5 FXS vs 5 control, `~ cell_type + condition` |
| Methylation | [GSE41273](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41273) | Peripheral blood, Illumina 450k | 9 case vs 53 control |

Both RNA-seq cohorts use **NCBI-reprocessed gene-level counts** (GRCh38.p13,
Entrez IDs), so quantification differences cannot account for non-replication
between them.

## Running (v5)

### Google Colab (recommended)
Open `colab/FXS_GRN_v5_Colab.ipynb` in Colab (Python runtime) and run the cells top to
bottom. It installs R 4.6 / Bioconductor 3.23 as binaries (r2u), fetches every input,
and writes data, results, checkpoints and logs to Google Drive. Long permutation steps
are resumable after a disconnect.

### Locally
```bash
export FXS_REPO=$PWD FXS_ROOT=$PWD      # FXS_ROOT = where outputs go
Rscript scripts/01_fetch_data.R         # downloads all inputs once
Rscript scripts/part1_steps1-8.R        # -> checkpoints/checkpoint_step8.RData
Rscript scripts/part2_steps9-16.R       # -> checkpoint_step16.RData
Rscript scripts/part3_steps17-19.R      # resumable nulls -> checkpoint_step19.RData
Rscript scripts/17b-run.R               # PANDA null, resumable
Rscript scripts/part4_steps20-22.R      # tissue generality, CollecTRI, summary
Rscript scripts/14_cytoscape_local.R    # with Cytoscape open: styled figures + .cys
```
Each part is a fresh R process (memory is released between them). Settings live in
`scripts/00_config.R` and can be overridden with environment variables. See
`CHANGES.md` for what changed from v4 and why.

**Stop at the FMR1 gate in Step 4.** If FMR1 is not significantly downregulated in FXS,
the sample labels are wrong and nothing downstream is meaningful.

## Pipeline

| Steps | Contents |
|---|---|
| 1–4 | Count loading, metadata-driven sample annotation, VST, DESeq2 |
| 5 | Differential methylation (limma), CpG→gene mapping, continuous promoter beta |
| 6–7 | TF list preparation, GENIE3 network inference |
| 8 | STRING PPI, JASPAR motif prior, methylation weighting, PANDA |
| 9–10 | igraph network construction, hub identification |
| 11–12 | LIONESS sample-specific networks, ALPACA differential modules |
| 13–14 | GO/KEGG enrichment, Cytoscape export |
| 15–16 | Validation cohort DE, cross-cohort concordance (Fisher, RRHO, effect-size correlation) |
| 17 | Permutation nulls: GENIE3, PANDA, ALPACA, LIONESS bootstrap |
| 18 | Weight-based centrality, target-restriction bias, module-local centrality |
| 19 | In silico validation: ChEA3, ENCODE ChIP-seq, TRRUST |

### Methylation integration

Promoter methylation enters the PANDA motif prior as a continuous weight
(`1.5 − promoter_beta`), computed from control-group beta values so the prior
encodes baseline promoter accessibility rather than the FXS contrast. An earlier
binary hyper/hypo scheme modulated only 0.12% of prior edges and was replaced.

Methylation comes from an independent cohort with no shared subjects. It is a
**population-level prior on regulatory accessibility, not matched per-subject
multi-omics.**

## Key results

**FMR1 silencing replicates across all three datasets:**

| Dataset | Result |
|---|---|
| GSE146339 | log2FC −10.47, padj 3.1 × 10⁻⁵¹ |
| GSE218676 | log2FC −8.91, padj 7.2 × 10⁻⁸ |
| GSE41273 | Hypermethylated, mean logFC 0.562 across 15 probes |

FMR1-AS1 is also downregulated (log2FC −5.08). Methylation is focal: FMR1 and
ASFMR1 account for the majority of significant DMPs.

**Known FXS regulators are locally central without being global hubs.** None of
the six appear in the GENIE3 genome-wide top 50, but 4 of 5 fall in the top
decile of weighted out-strength within their own ALPACA module (median
percentile 0.953).

**In silico validation:** ENCODE ChIP-seq self-recovery places 6/18 TF-network
pairs at rank ≤10 and 9/18 at ≤25 (REST and NRF1 both rank 1). TRRUST curated
targets are enriched for 3/20 pairs at FDR < 0.05, led by SP1 (4.3× enrichment,
p = 1.8 × 10⁻⁷).

## Negative results

These are reported deliberately. Four independent tests agree that network-level
inference is underpowered at this sample size.

| Test | Result |
|---|---|
| Cross-cohort effect-size correlation | Pearson r = 0.034, sign concordance 50.6% |
| DEG overlap | 1 gene (FMR1), but 94.7× enriched over chance (p = 0.011) |
| GENIE3 permutation null | 46,796 edges nominal p < 0.05, **0 at FDR < 0.05** |
| PANDA randomized-prior null | Null q95 = 3.49 vs observed 90th pct 3.51 — indistinguishable |
| ALPACA permutation null | 484 nodes nominal, **0 at FDR**; best module FDR = 0.19 |
| LIONESS bootstrap | 49.9% of edges stronger in FXS, 95% CI [46.7%, 52.3%] |

The PANDA result specifically indicates that edge scores largely reproduce the
motif/PPI prior rather than adding expression-driven signal.

## Limitations

1. Cohorts differ in tissue (cerebral organoid vs NPC / B-lymphocyte / iPSC).
   Genes replicating across both are robust to tissue context; failures to
   replicate may be tissue-specific rather than false positives.
2. Methylation is from a third tissue (blood) and a third cohort.
3. GENIE3 ran on 1,340 genes of which 1,330 are TFs — only 13 DEGs were
   available to seed the target set. It is effectively a TF–TF network, and
   only ~0.3% of each known TF's motif regulon was eligible to be discovered
   (see `TF_target_eligibility_bias.csv`).
4. LIONESS has no covariate mechanism, so developmental timepoint is unadjusted
   in the sample-specific networks.
5. Validation cohort: NCBI reprocessed 20 of 26 submitted samples (FXS_426_NPC
   absent). Two donors each contribute two cell types, so effective n is below
   the line count.
6. In silico validation shows inferred edges recover known binding above chance.
   It does not establish causality in FXS and is not equivalent to experimental
   confirmation.

**Note on v4 GENIE3 outputs:** in v4 the GENIE3 regulator and target columns were
swapped (GENIE3 normalises per target, hence `Sum_Weight = 1` per "TF"). Fixed in v5;
see `CHANGES.md`.

## Environment

R 4.6, Bioconductor 3.23 (exact versions written to `logs/package_versions.csv` and
`logs/sessionInfo_*.txt` on every run). Core packages: DESeq2, limma, minfi, GENIE3, pandaR,
netZooR, igraph, clusterProfiler, JASPAR2020, motifmatchr,
BSgenome.Hsapiens.UCSC.hg38, RCy3. Run `sessionInfo()` for exact versions.

netZooR/pandaR in this build contain `class(x) == "ExpressionSet"` comparisons
that fail under R ≥ 4.0; Step 11 patches these at the namespace level before
calling LIONESS.
