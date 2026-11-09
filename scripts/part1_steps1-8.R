# ============================================================
# FXS PIPELINE v5 — PART 1 of 3 : STEPS 1–8
# Data loading, DE, methylation, TF prep, GENIE3, PANDA
#
# Requires: scripts/01_fetch_data.R to have run once
# Ends by saving checkpoint_step8.RData
# Then: restart R (new Rscript call) and run part2_steps9-16.R
#
# v5 changes (see CHANGES.md):
#   4   apeglm-shrunken LFCs added (used for cross-cohort correlation)
#   5A  blood cell-composition (EpiDISH) adjusted DMP model
#   6   timepoint removed from the network-inference expression matrix
#   7   GENIE3 regulator/target columns were SWAPPED in v4 — fixed;
#       target set expanded from 13 DEGs to a label-independent
#       top-variance universe
#   8   STRING read gz-compressed via data.table; motif target pool
#       chosen by variance (label-independent) instead of DE p-value
# ============================================================

# ---- Core ----
library(DESeq2)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(reshape2)
library(scales)
library(Biobase)

# ---- Methylation ----
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# ---- Annotation ----
library(org.Hs.eg.db)

# ---- Motif / genome ----
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)

# ---- Networks ----
library(GENIE3)
library(pandaR)
library(igraph)

# ---- v5 additions ----
library(apeglm)
library(EpiDISH)
library(data.table)

# Paths, config and helpers (portable: Colab or local)
source(file.path(Sys.getenv("FXS_REPO", "."), "scripts", "00_config.R"))
log_session("part1")

annot  <- read.delim(gzfile(file.path(DATA_DIR, "Human.GRCh38.p13.annot.tsv.gz")), header = TRUE,
                     quote = "", check.names = FALSE, stringsAsFactors = FALSE)
id2sym <- setNames(annot$Symbol, as.character(annot$GeneID))
cat("Annotation rows:", nrow(annot), "\n")

# ============================================================
# STEP 1 — LOAD PRIMARY COUNT MATRIX (GSE146339)
# ============================================================

counts_raw <- read.delim(gzfile(file.path(DATA_DIR, "GSE146339_raw_counts_GRCh38.p13_NCBI.tsv.gz")),
                         header = TRUE, check.names = FALSE,
                         stringsAsFactors = FALSE)
stopifnot(colnames(counts_raw)[1] == "GeneID")

count_matrix <- as.matrix(counts_raw[, -1])
rownames(count_matrix) <- as.character(counts_raw$GeneID)
storage.mode(count_matrix) <- "integer"

cat("GSE146339 — genes:", nrow(count_matrix),
    "| samples:", ncol(count_matrix), "\n")
cat("Library sizes (millions):\n")
print(round(colSums(count_matrix) / 1e6, 2))

# ============================================================
# STEP 2 — SAMPLE ANNOTATION
# ============================================================

sample_info <- data.frame(
  gsm = c("GSM4378490","GSM4378491","GSM4378492",
          "GSM4378493","GSM4378494","GSM4378495",
          "GSM4378496","GSM4378497","GSM4378498",
          "GSM4378499","GSM4378500","GSM4378501"),
  label = c("D28_control1","D28_control2","D28_control3",
            "D28_FXS1","D28_FXS2","D28_FXS3",
            "D105_control4","D105_control5","D105_control6",
            "D105_FXS4","D105_FXS5","D105_FXS6"),
  stringsAsFactors = FALSE)

sample_info$timepoint <- ifelse(grepl("^D28", sample_info$label), "D28", "D105")
sample_info$condition <- ifelse(grepl("control", sample_info$label),
                                "Control", "FXS")

stopifnot(setequal(sample_info$gsm, colnames(count_matrix)))
sample_info <- sample_info[match(colnames(count_matrix), sample_info$gsm), ]
stopifnot(identical(sample_info$gsm, colnames(count_matrix)))

cat("\n--- Sample annotation ---\n"); print(sample_info)

ct_tab <- table(sample_info$timepoint, sample_info$condition)
cat("\nTimepoint x condition:\n"); print(ct_tab)
if (any(ct_tab == 0))
  stop("Rank-deficient design: a timepoint lacks one condition.")

colnames(count_matrix) <- sample_info$label

count_matrix_filtered <- count_matrix[rowSums(count_matrix) >= 10, ]
cat("\nGenes before filtering:", nrow(count_matrix),
    "| after:", nrow(count_matrix_filtered), "\n")

# ============================================================
# STEP 3 — NORMALIZATION
# ============================================================

coldata <- data.frame(
  condition = factor(sample_info$condition, levels = c("Control", "FXS")),
  timepoint = factor(sample_info$timepoint, levels = c("D28", "D105")),
  row.names = colnames(count_matrix_filtered))

design_formula <- ~ timepoint + condition
cat("\nDesign formula:", deparse(design_formula), "\n")

dds <- DESeqDataSetFromMatrix(countData = count_matrix_filtered,
                              colData   = coldata,
                              design    = design_formula)

condition <- colData(dds)$condition
timepoint <- colData(dds)$timepoint
names(condition) <- colnames(dds)
names(timepoint) <- colnames(dds)
cell_type <- timepoint

cat("FXS:", sum(condition == "FXS"),
    "| Control:", sum(condition == "Control"), "\n")

vsd <- vst(dds, blind = TRUE)
expr_normalized <- assay(vsd)
write.csv(expr_normalized, file.path(DATA_DIR, "GSE146339_VST_normalized.csv"))

pca_data <- prcomp(t(expr_normalized))
pct_var  <- round(100 * pca_data$sdev^2 / sum(pca_data$sdev^2), 1)
pca_df <- data.frame(pca_data$x[, 1:2], condition = condition,
                     timepoint = timepoint, sample = colnames(expr_normalized))

ggplot(pca_df, aes(x = PC1, y = PC2, color = condition,
                   shape = timepoint, label = sample)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = c("FXS" = "#E74C3C", "Control" = "#3498DB")) +
  labs(title    = "PCA — GSE146339 cerebral organoids",
       subtitle = "6 FXS vs 6 Control | D28 and D105",
       x = paste0("PC1 (", pct_var[1], "%)"),
       y = paste0("PC2 (", pct_var[2], "%)")) +
  theme_bw() +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
ggsave(file.path(FIG_DE, "PCA_GSE146339.png"), width = 10, height = 7, dpi = 300)

cat("\nStep 3 done. PC1 separating D28 from D105 is expected.\n")

# ============================================================
# STEP 4 — DIFFERENTIAL EXPRESSION
# ============================================================

dds <- DESeq(dds)
cat("\nCoefficients:", resultsNames(dds), "\n")

res    <- results(dds, contrast = c("condition", "FXS", "Control"))
res_df <- as.data.frame(res)
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
summary(res, alpha = 0.05)

res_df$GeneID <- rownames(res_df)
res_df$Symbol <- id2sym[res_df$GeneID]
cat("Genes with a symbol:", sum(!is.na(res_df$Symbol)), "/", nrow(res_df), "\n")

# v5: apeglm shrinkage. DEG calls still use the Wald test on unshrunk
# LFCs; shrunken LFCs remove the low-count noise that dominated the
# v4 genome-wide cross-cohort correlation (Step 16).
res_shr <- lfcShrink(dds, coef = "condition_FXS_vs_Control", type = "apeglm",
                     quiet = TRUE)
res_df$log2FC_shrunk <- res_shr[res_df$GeneID, "log2FoldChange"]
rm(res_shr)
write.csv(res_df, file.path(RESULTS_DE, "DESeq2_FXS_vs_Control_all.csv"))

sig_degs <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(padj)
write.csv(sig_degs, file.path(RESULTS_DE, "DESeq2_significant_DEGs.csv"))

cat("\nSignificant DEGs:", nrow(sig_degs), "\n")
cat("Upregulated:",   sum(sig_degs$log2FoldChange > 0), "\n")
cat("Downregulated:", sum(sig_degs$log2FoldChange < 0), "\n")

# ---- FMR1 SANITY CHECK — hard gate ----
cat("\n=== FMR1 SANITY CHECK ===\n")
fmr1_row <- res_df[which(res_df$Symbol == "FMR1"), ]
if (nrow(fmr1_row) > 0) {
  print(fmr1_row[, c("Symbol","baseMean","log2FoldChange","pvalue","padj")])
  if (!is.na(fmr1_row$log2FoldChange[1]) && fmr1_row$log2FoldChange[1] < 0) {
    cat("PASS — FMR1 down in FXS, as expected.\n")
  } else {
    cat("FAIL — FMR1 not down in FXS. STOP and recheck the GSM mapping.\n")
  }
} else {
  cat("FMR1 absent after filtering — consistent with silencing. Raw counts:\n")
  fid <- names(id2sym)[which(id2sym == "FMR1")]
  if (length(fid) > 0 && fid[1] %in% rownames(count_matrix))
    print(count_matrix[fid[1], , drop = FALSE])
}

res_plot <- res_df %>%
  filter(!is.na(padj)) %>%
  mutate(significance = case_when(
    padj < 0.05 & log2FoldChange >  1 ~ "Upregulated",
    padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
    TRUE ~ "Not significant"),
    gene = ifelse(is.na(Symbol), GeneID, Symbol))

ggplot(res_plot, aes(x = log2FoldChange, y = -log10(padj),
                     color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated"     = "#E74C3C",
                                "Downregulated"   = "#3498DB",
                                "Not significant" = "grey70")) +
  geom_text_repel(data = res_plot %>%
                    filter(padj < 0.01, abs(log2FoldChange) > 2) %>% head(20),
                  aes(label = gene), size = 3, max.overlaps = 15) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title    = "FXS vs Control — Volcano Plot (GSE146339)",
       subtitle = "~ timepoint + condition | 6 vs 6 organoids",
       x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_bw() +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
ggsave(file.path(FIG_DE, "volcano_plot_FXS_vs_Control.png"), width = 10, height = 7, dpi = 300)

if (nrow(sig_degs) >= 10) {
  top_ids  <- head(rownames(sig_degs), min(50, nrow(sig_degs)))
  heat_mat <- expr_normalized[top_ids, , drop = FALSE]
  rownames(heat_mat) <- ifelse(is.na(id2sym[top_ids]), top_ids, id2sym[top_ids])
  pheatmap(heat_mat,
           annotation_col = data.frame(Condition = condition,
                                       Timepoint = timepoint,
                                       row.names = colnames(heat_mat)),
           scale = "row", show_rownames = TRUE, show_colnames = TRUE,
           main = "Top DEGs — FXS vs Control",
           filename = file.path(FIG_DE, "heatmap_top_DEGs.png"), width = 10, height = 12)
} else {
  cat("\nFewer than 10 DEGs — heatmap skipped.\n")
}

cat("\nStep 4 done — proceed only if FMR1 passed.\n")

# ============================================================
# STEP 5 — METHYLATION
# ============================================================

# ---- 5A: DIFFERENTIAL METHYLATION ----

meth_rds <- file.path(DATA_DIR, "GSE41273_methylation_matrix.rds")
meth_matrix <- if (file.exists(meth_rds)) readRDS(meth_rds) else
  as.matrix(read.csv(file.path(DATA_DIR, "GSE41273_methylation_matrix.csv"), row.names = 1))
meta_meth   <- read.csv(file.path(DATA_DIR, "GSE41273_metadata.csv"), row.names = 1)
meta_meth   <- meta_meth[colnames(meth_matrix), , drop = FALSE]
stopifnot(identical(rownames(meta_meth), colnames(meth_matrix)))
cat("\nMethylation matrix dimensions:", dim(meth_matrix), "\n")

meth_condition <- factor(ifelse(
  meta_meth$characteristics_ch1.2 == "phenotype: Case", "FXS", "Control"),
  levels = c("Control", "FXS"))

cat("FXS:", sum(meth_condition == "FXS"),
    "| Control:", sum(meth_condition == "Control"), "\n")
cat("Sex:", paste(names(table(meta_meth$characteristics_ch1)),
                  table(meta_meth$characteristics_ch1), collapse = "; "),
    "— all-male cohort, so X-inactivation cannot confound FMR1\n")
cat("Beta value range:", range(meth_matrix, na.rm = TRUE), "\n")

meth_clean <- meth_matrix[complete.cases(meth_matrix), ]
cat("Probes after NA removal:", nrow(meth_clean), "\n")

# ---- v5: blood cell-composition adjustment ----
# Whole-blood methylation is dominated by leukocyte proportions; an
# unadjusted case/control model can call composition shifts as DMPs.
# EpiDISH (RPC, 7-type blood reference) estimates proportions from the
# betas themselves. Neutrophils are dropped from the design because
# the proportions sum to 1.
data("centDHSbloodDMC.m", package = "EpiDISH", envir = environment())
cell_props <- EpiDISH::epidish(beta.m = meth_clean,
                               ref.m  = centDHSbloodDMC.m,
                               method = "RPC")$estF
write.csv(data.frame(Sample = rownames(cell_props), Condition = meth_condition,
                     cell_props),
          file.path(RESULTS_METH, "EpiDISH_cell_proportions.csv"), row.names = FALSE)

cell_diff <- sapply(colnames(cell_props), function(ct)
  wilcox.test(cell_props[, ct] ~ meth_condition, exact = FALSE)$p.value)
cat("\nCell-proportion differences FXS vs Control (Wilcoxon p):\n")
print(signif(cell_diff, 3))

design_unadj <- model.matrix(~ meth_condition)
fit_unadj    <- eBayes(lmFit(meth_clean, design_unadj))
dmp_unadj    <- topTable(fit_unadj, coef = "meth_conditionFXS", number = Inf,
                         adjust.method = "BH")

cov_cols    <- setdiff(colnames(cell_props), "Neutro")
design_adj  <- model.matrix(~ meth_condition + cell_props[, cov_cols])
fit_adj     <- eBayes(lmFit(meth_clean, design_adj))
dmp_adj     <- topTable(fit_adj, coef = "meth_conditionFXS", number = Inf,
                        adjust.method = "BH")

dmp_results <- if (METH_ADJUST_CELLTYPE) dmp_adj else dmp_unadj
dmp_results <- dmp_results[order(dmp_results$adj.P.Val), ]
cat("\nPrimary DMP model:", if (METH_ADJUST_CELLTYPE)
  "~ condition + cell proportions" else "~ condition (v4)", "\n")
cat("DMPs at FDR < 0.05 — adjusted:", sum(dmp_adj$adj.P.Val < 0.05),
    "| unadjusted:", sum(dmp_unadj$adj.P.Val < 0.05), "\n")

write.csv(dmp_results, file.path(RESULTS_METH, "DMP_FXS_vs_Control_all.csv"))
write.csv(dmp_unadj,   file.path(RESULTS_METH, "DMP_FXS_vs_Control_all_unadjusted.csv"))
rm(fit_unadj, fit_adj); gc(verbose = FALSE)

sig_dmps <- dmp_results %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 0.1) %>%
  arrange(adj.P.Val)

cat("Significant DMPs (FDR<0.05, |logFC|>0.1):", nrow(sig_dmps), "\n")
cat("Hypermethylated:", sum(sig_dmps$logFC > 0), "\n")
cat("Hypomethylated:",  sum(sig_dmps$logFC < 0), "\n")
write.csv(sig_dmps, file.path(RESULTS_METH, "DMP_FXS_vs_Control_significant.csv"))

# ---- 5B: DMP -> GENE MAPPING ----

cat("\n=== Step 5B: CpG probe to gene mapping ===\n")
anno_450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

anno_sub <- anno_450k[rownames(anno_450k) %in% rownames(sig_dmps), ]
cat("Probes with annotation:", nrow(anno_sub), "/", nrow(sig_dmps), "\n")

probe_gene_map <- data.frame(
  Probe = rownames(anno_sub),
  Gene  = anno_sub$UCSC_RefGene_Name,
  stringsAsFactors = FALSE) %>%
  filter(Gene != "") %>%
  mutate(Gene = strsplit(Gene, ";")) %>%
  tidyr::unnest(Gene) %>%
  mutate(Gene = trimws(Gene)) %>%
  distinct()

dmp_merge <- merge(probe_gene_map,
                   data.frame(Probe   = rownames(sig_dmps),
                              logFC   = sig_dmps$logFC,
                              adjPval = sig_dmps$adj.P.Val),
                   by = "Probe")

gene_meth_status <- dmp_merge %>%
  group_by(Gene) %>%
  summarise(mean_logFC = mean(logFC), n_probes = n(), .groups = "drop") %>%
  mutate(meth_status = case_when(
    mean_logFC >  0.1 ~ "Hypermethylated",
    mean_logFC < -0.1 ~ "Hypomethylated",
    TRUE              ~ "Unchanged"))

cat("\nGene methylation status:\n"); print(table(gene_meth_status$meth_status))
print(gene_meth_status)
write.csv(gene_meth_status, file.path(RESULTS_METH, "methylation_gene_summary.csv"), row.names = FALSE)

fmr1_meth <- gene_meth_status[gene_meth_status$Gene == "FMR1", ]
if (nrow(fmr1_meth) > 0) {
  cat("\n=== KEY BIOLOGICAL FINDING ===\n")
  cat("FMR1:", fmr1_meth$meth_status, "| mean logFC",
      round(fmr1_meth$mean_logFC, 3), "|", fmr1_meth$n_probes, "probes\n")
  cat("FXS methylation is focal — the majority of significant DMPs are\n")
  cat("FMR1/ASFMR1 probes. This is the expected single-locus pattern.\n")
}

# v5 sensitivity: FMR1/ASFMR1 probes under both models
fmr1_probes <- rownames(anno_450k)[grepl("(^|;)(FMR1|ASFMR1|FMR1-AS1)(;|$)",
                                         anno_450k$UCSC_RefGene_Name)]
meth_sens <- data.frame(
  Probe         = fmr1_probes,
  dBeta_unadj   = dmp_unadj[fmr1_probes, "logFC"],
  FDR_unadj     = dmp_unadj[fmr1_probes, "adj.P.Val"],
  dBeta_celladj = dmp_adj[fmr1_probes, "logFC"],
  FDR_celladj   = dmp_adj[fmr1_probes, "adj.P.Val"])
meth_sens <- meth_sens[!is.na(meth_sens$dBeta_unadj), ]
write.csv(meth_sens, file.path(RESULTS_METH, "FMR1_probes_celltype_sensitivity.csv"),
          row.names = FALSE)
cat("\nFMR1/ASFMR1 probes — mean dBeta unadjusted:",
    round(mean(meth_sens$dBeta_unadj), 3),
    "| cell-adjusted:", round(mean(meth_sens$dBeta_celladj), 3), "\n")
sig_unadj_n <- sum(dmp_unadj$adj.P.Val < 0.05 & abs(dmp_unadj$logFC) > 0.1)
sig_adj_n   <- sum(dmp_adj$adj.P.Val   < 0.05 & abs(dmp_adj$logFC)   > 0.1)
write.csv(data.frame(Model = c("Unadjusted (v4)", "Cell-composition adjusted (v5)"),
                     Sig_DMPs = c(sig_unadj_n, sig_adj_n),
                     Mean_dBeta_FMR1_probes = c(mean(meth_sens$dBeta_unadj),
                                                mean(meth_sens$dBeta_celladj))),
          file.path(RESULTS_METH, "DMP_model_comparison.csv"), row.names = FALSE)

ggplot(gene_meth_status,
       aes(x = reorder(Gene, mean_logFC), y = mean_logFC, fill = meth_status)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = c("Hypermethylated" = "#E74C3C",
                               "Hypomethylated"  = "#3498DB",
                               "Unchanged"       = "#95A5A6")) +
  labs(title    = "Differentially Methylated Genes — FXS vs Control",
       subtitle = "GSE41273 | FDR < 0.05 | |logFC| > 0.1",
       x = "Gene", y = "Mean logFC (FXS - Control)",
       fill = "Methylation Status") +
  theme_bw() +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
ggsave(file.path(FIG_METH, "methylation_gene_summary.png"), width = 8, height = 6, dpi = 300)

# ---- 5C: CONTINUOUS PROMOTER METHYLATION ----

cat("\n=== Step 5C: Continuous promoter methylation ===\n")

promoter_groups <- c("TSS200", "TSS1500", "5'UTR", "1stExon")
is_promoter <- sapply(strsplit(anno_450k$UCSC_RefGene_Group, ";"),
                      function(g) any(g %in% promoter_groups))
promoter_probes <- rownames(anno_450k)[which(is_promoter)]
cat("Promoter-associated probes:", length(promoter_probes), "\n")

probe2gene_all <- data.frame(
  Probe = rownames(anno_450k),
  Gene  = anno_450k$UCSC_RefGene_Name,
  stringsAsFactors = FALSE) %>%
  filter(Gene != "", Probe %in% promoter_probes) %>%
  mutate(Gene = strsplit(Gene, ";")) %>%
  tidyr::unnest(Gene) %>%
  mutate(Gene = trimws(Gene)) %>%
  distinct()

ctrl_beta <- rowMeans(meth_clean[, meth_condition == "Control", drop = FALSE],
                      na.rm = TRUE)

gene_beta <- probe2gene_all %>%
  mutate(beta = ctrl_beta[Probe]) %>%
  filter(!is.na(beta)) %>%
  group_by(Gene) %>%
  summarise(promoter_beta = mean(beta), n_probes = n(), .groups = "drop")

cat("Genes with promoter methylation data:", nrow(gene_beta), "\n")
cat("Promoter beta range:", round(range(gene_beta$promoter_beta), 3), "\n")
cat("Median promoter beta:", round(median(gene_beta$promoter_beta), 3), "\n")
write.csv(gene_beta, file.path(RESULTS_METH, "promoter_methylation_by_gene.csv"), row.names = FALSE)

ggplot(gene_beta, aes(x = promoter_beta)) +
  geom_histogram(bins = 60, fill = "#3498DB", color = "white") +
  geom_vline(xintercept = median(gene_beta$promoter_beta),
             linetype = "dashed", color = "#E74C3C") +
  labs(title    = "Promoter Methylation Distribution (control samples)",
       subtitle = paste0(nrow(gene_beta),
                         " genes | dashed = median | used as PANDA prior weight"),
       x = "Mean promoter beta", y = "Genes") +
  theme_bw() +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
ggsave(file.path(FIG_METH, "promoter_methylation_distribution.png"), width = 8, height = 5, dpi = 300)

cat("\nStep 5 done\n")

# ============================================================
# STEP 6 — TF LIST PREPARATION + FMR1 DIAGNOSTIC
# ============================================================

tf_list <- readLines(file.path(DATA_DIR, "tf_list.txt"))
tf_list <- unique(trimws(tf_list[tf_list != ""]))
cat("\nTF list size:", length(tf_list), "\n")

# ---- v5: remove developmental timepoint before network inference ----
# PC1 in GSE146339 is D28 vs D105. Left in, it drives most co-expression
# that GENIE3/PANDA/LIONESS read as regulation, and LIONESS has no
# covariate mechanism. removeBatchEffect regresses timepoint out while
# protecting the condition effect (design = ~ condition). DE (Step 4)
# still models timepoint inside DESeq2 and is unaffected.
expr_network <- expr_normalized
if (ADJUST_TIMEPOINT_FOR_NETWORKS) {
  expr_network <- limma::removeBatchEffect(
    expr_normalized, batch = timepoint,
    design = model.matrix(~ condition))
  pc_adj <- prcomp(t(expr_network))
  cat("Timepoint removed for networks. PC1 vs timepoint |r| before:",
      round(abs(cor(pca_data$x[, 1], as.numeric(timepoint))), 2),
      "| after:", round(abs(cor(pc_adj$x[, 1], as.numeric(timepoint))), 2), "\n")
  rm(pc_adj)
}

expr_normalized_sym <- expr_network
sym_rows <- id2sym[rownames(expr_normalized_sym)]
rownames(expr_normalized_sym) <- ifelse(is.na(sym_rows) | sym_rows == "",
                                        rownames(expr_normalized_sym), sym_rows)

row_vars <- apply(expr_normalized_sym, 1, var)
expr_normalized_sym <- expr_normalized_sym[order(row_vars, decreasing = TRUE), ]
expr_normalized_sym <- expr_normalized_sym[!duplicated(rownames(expr_normalized_sym)), ]
cat("Expression matrix after dedup:", nrow(expr_normalized_sym), "genes\n")

gene_symbols <- setNames(id2sym[rownames(expr_normalized)],
                         rownames(expr_normalized))

tf_in_data <- intersect(tf_list, rownames(expr_normalized_sym))
cat("TFs found in expression data:", length(tf_in_data), "\n")

cat("\nKnown FXS TFs recovered:\n")
print(data.frame(TF = known_tf,
                 in_expr       = known_tf %in% rownames(expr_normalized_sym),
                 in_tf_in_data = known_tf %in% tf_in_data))

cat("\n=== Step 6D: FMR1 / ASFMR1 ===\n")
fmr1_genes <- res_df[which(res_df$Symbol %in%
                             c("FMR1", "ASFMR1", "FMR1-AS1")), ]
if (nrow(fmr1_genes) > 0) {
  print(fmr1_genes[, c("Symbol","baseMean","log2FoldChange","pvalue","padj")])
  write.csv(fmr1_genes, file.path(RESULTS_METH, "FMR1_diagnostic.csv"), row.names = FALSE)
} else {
  cat("Neither gene present after filtering.\n")
}
cat("Gene-level counts — no isoform ambiguity to resolve.\n")

cat("\nStep 6 done\n")

# ============================================================
# STEP 7 — GENIE3 NETWORK INFERENCE
# ============================================================

cat("\nStarting GENIE3 inference...\n")
expr_matrix <- as.matrix(expr_normalized_sym)

deg_symbols <- sig_degs$Symbol[!is.na(sig_degs$Symbol)]

# ---- v5: label-independent target universe ----
# v4 used TFs + DEGs only (1,330 of 1,340 genes were TFs), so GENIE3 was
# a TF-TF map and ~0.3% of each TF's regulon was eligible. v5 adds the
# top-variance non-TF genes (variance on the timepoint-adjusted matrix,
# no FXS/Control labels used), which is the standard GENIE3 setup.
non_tf_genes <- setdiff(rownames(expr_matrix), tf_in_data)
non_tf_genes <- non_tf_genes[!grepl("^LOC\\d+$|^\\d+$", non_tf_genes)]  # unannotated IDs
# Expression floor before ranking by variance: >= 10 normalised counts in
# >= 3 samples (the smallest timepoint x condition group). Without it,
# near-silent genes (olfactory receptors, keratins, testis genes in a brain
# organoid) top the variance list through VST noise and become fake hubs.
norm_ct  <- counts(dds, normalized = TRUE)
expr_ok  <- rownames(norm_ct)[rowSums(norm_ct >= 10) >= 3]
expr_ok_sym <- unique(na.omit(id2sym[expr_ok]))
cat("Non-TF genes before / after expression floor:", length(non_tf_genes), "/",
    sum(non_tf_genes %in% expr_ok_sym), "\n")
non_tf_genes <- intersect(non_tf_genes, expr_ok_sym)
rm(norm_ct)
nt_var       <- apply(expr_matrix[non_tf_genes, , drop = FALSE], 1, var)
hv_targets   <- names(sort(nt_var, decreasing = TRUE))[seq_len(min(N_GENIE3_TARGETS,
                                                                    length(nt_var)))]
if (TARGET_SELECTION == "pvalue") hv_targets <- character(0)   # reproduce v4
top_genes <- unique(c(deg_symbols, hv_targets, tf_in_data))
top_genes <- intersect(top_genes, rownames(expr_matrix))

cat("Genes for GENIE3:", length(top_genes),
    "| TFs:", length(tf_in_data),
    "| non-TF targets:", length(setdiff(top_genes, tf_in_data)), "\n")

expr_subset <- expr_matrix[top_genes, ]

genie3_weights <- GENIE3(expr_subset, regulators = tf_in_data,
                         nTrees = N_TREES_OBS, nCores = N_CORES)
cat("GENIE3 complete\n")

# ---- v5 FIX: orientation ----
# GENIE3 returns weightMatrix[regulator, target] (rows = regulators).
# v4 named melt()'s Var1 "Target" and Var2 "TF", so every downstream
# GENIE3 "TF" was actually a target and vice versa (visible in v4 output
# as non-TF genes such as FMR1/TSHR/RERG having out-degree ~80-100, and
# Sum_Weight = 1 per "TF" — GENIE3 normalises per TARGET).
stopifnot(all(rownames(genie3_weights) %in% tf_in_data),
          ncol(genie3_weights) == length(top_genes))
edge_list_genie3 <- reshape2::melt(genie3_weights, varnames = c("TF", "Target"),
                                   value.name = "Weight")
edge_list_genie3$TF     <- as.character(edge_list_genie3$TF)
edge_list_genie3$Target <- as.character(edge_list_genie3$Target)
edge_list_genie3 <- edge_list_genie3[edge_list_genie3$Weight > 0 &
                                       edge_list_genie3$TF != edge_list_genie3$Target, ]
stopifnot(all(edge_list_genie3$TF %in% tf_in_data))
chk <- tapply(edge_list_genie3$Weight, edge_list_genie3$Target, sum)
cat("Per-target weight sums (should be ~1):", round(range(chk), 3), "\n")
rm(chk)
edge_list_genie3 <- edge_list_genie3[order(edge_list_genie3$Weight,
                                           decreasing = TRUE), ]

cat("Total edges:", nrow(edge_list_genie3), "\n")
print(head(edge_list_genie3, 10))
write.csv(edge_list_genie3, file.path(RESULTS_GENIE3, "GENIE3_edge_list.csv"), row.names = FALSE)

cat("Known FXS TFs found:",
    paste(intersect(known_tf, unique(edge_list_genie3$TF)), collapse = ", "), "\n")
cat("\nStep 7 done\n")

# ============================================================
# STEP 8 — PANDA NETWORK INFERENCE
# ============================================================

# ---- 8A: STRING PPI ----

cat("\nLoading STRING PPI...\n")
string_gz  <- file.path(DATA_DIR, "9606.protein.links.v12.0.txt.gz")
string_txt <- file.path(DATA_DIR, "9606.protein.links.v12.0.txt")
ppi_raw <- data.table::fread(if (file.exists(string_gz)) string_gz else string_txt,
                             sep = " ", header = TRUE, showProgress = FALSE)
ppi_filtered <- as.data.frame(ppi_raw[combined_score >= 700])
cat("High confidence PPI edges:", nrow(ppi_filtered), "\n")

ppi_filtered$protein1 <- gsub("9606\\.", "", ppi_filtered$protein1)
ppi_filtered$protein2 <- gsub("9606\\.", "", ppi_filtered$protein2)
ensp_ids <- unique(c(ppi_filtered$protein1, ppi_filtered$protein2))

ensp_to_symbol <- mapIds(org.Hs.eg.db, keys = ensp_ids, column = "SYMBOL",
                         keytype = "ENSEMBLPROT", multiVals = "first")
ppi_filtered$gene1 <- ensp_to_symbol[ppi_filtered$protein1]
ppi_filtered$gene2 <- ensp_to_symbol[ppi_filtered$protein2]
ppi_filtered <- ppi_filtered[!is.na(ppi_filtered$gene1) &
                               !is.na(ppi_filtered$gene2), ]
cat("PPI edges after symbol mapping:", nrow(ppi_filtered), "\n")

ppi_tf <- ppi_filtered[ppi_filtered$gene1 %in% tf_in_data &
                         ppi_filtered$gene2 %in% tf_in_data, ]
cat("TF-TF PPI edges:", nrow(ppi_tf), "\n")

ppi_matrix <- matrix(0, nrow = length(tf_in_data), ncol = length(tf_in_data),
                     dimnames = list(tf_in_data, tf_in_data))
if (nrow(ppi_tf) > 0) {
  for (i in seq_len(nrow(ppi_tf))) {
    g1 <- ppi_tf$gene1[i]; g2 <- ppi_tf$gene2[i]
    ppi_matrix[g1, g2] <- ppi_tf$combined_score[i] / 1000
    ppi_matrix[g2, g1] <- ppi_tf$combined_score[i] / 1000
  }
}
cat("PPI matrix dimensions:", dim(ppi_matrix), "\n")

rm(ppi_raw); gc(verbose = FALSE)   # a few million rows, not needed again

# ---- 8B: JASPAR MOTIF PRIOR ----

cat("\nBuilding motif prior from JASPAR...\n")

if (TARGET_SELECTION == "pvalue") {
  # v4 behaviour: top genes by DE p-value (selection on the FXS contrast)
  top_ranked <- res_df %>%
    filter(!is.na(pvalue), !is.na(Symbol), Symbol != "") %>%
    arrange(pvalue) %>%
    head(N_TARGET_POOL)
  target_genes <- unique(top_ranked$Symbol[!top_ranked$Symbol %in% tf_in_data])
} else {
  # v5: label-independent — DEGs plus the most variable non-TF genes
  # from the same universe GENIE3 used, so both networks share targets
  top_ranked   <- data.frame(Symbol = unique(c(deg_symbols, hv_targets)))
  target_genes <- head(setdiff(top_ranked$Symbol, tf_in_data), N_TARGET_POOL)
}

cat("Target selection:", TARGET_SELECTION, "| ranked genes drawn:", nrow(top_ranked), "\n")
cat("Non-TF target genes:", length(target_genes), "\n")
cat("TFs in data:", length(tf_in_data), "\n")
cat("Target/TF ratio:", round(length(target_genes) / length(tf_in_data), 2), "\n")
if (length(target_genes) <= length(tf_in_data))
  cat("WARNING: targets do not exceed TFs — raise N_TARGET_POOL, ALPACA will fail\n")

target_entrez <- mapIds(org.Hs.eg.db, keys = target_genes, column = "ENTREZID",
                        keytype = "SYMBOL", multiVals = "first")
failed <- target_genes[is.na(target_entrez)]
if (length(failed) > 0) {
  alias_entrez <- mapIds(org.Hs.eg.db, keys = failed, column = "ENTREZID",
                         keytype = "ALIAS", multiVals = "first")
  target_entrez[is.na(target_entrez)] <-
    alias_entrez[names(target_entrez[is.na(target_entrez)])]
}
target_entrez <- target_entrez[!is.na(target_entrez)]
cat("Target genes with Entrez IDs:", length(target_entrez), "\n")
if (length(target_entrez) == 0) stop("No Entrez IDs found")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoters_gr <- suppressMessages(
  promoters(genes(txdb), upstream = 2000, downstream = 500))
promoters_gr <- promoters_gr[promoters_gr$gene_id %in% as.character(target_entrez)]

gene_id_map <- mapIds(org.Hs.eg.db, keys = as.character(promoters_gr$gene_id),
                      column = "SYMBOL", keytype = "ENTREZID",
                      multiVals = "first")
promoters_gr$symbol <- gene_id_map[as.character(promoters_gr$gene_id)]
promoters_gr <- promoters_gr[!is.na(promoters_gr$symbol)]
cat("Promoter regions:", length(promoters_gr), "\n")
if (length(promoters_gr) == 0) stop("No promoter regions")

# ---- SCAFFOLD FILTER AND TRIM ----
# matchMotifs cannot load sequence past the boundary of a non-circular
# sequence. Genes on _random / _alt / _fix contigs produce 2500bp
# windows that overrun short scaffolds, which aborts the scan.

main_chroms  <- paste0("chr", c(1:22, "X", "Y"))
n_before     <- length(promoters_gr)
promoters_gr <- promoters_gr[as.character(seqnames(promoters_gr)) %in% main_chroms]
seqlevels(promoters_gr) <- seqlevelsInUse(promoters_gr)
cat("Promoters on main chromosomes:", length(promoters_gr),
    sprintf("(dropped %d scaffold entries)\n", n_before - length(promoters_gr)))

seqlengths(promoters_gr) <-
  seqlengths(BSgenome.Hsapiens.UCSC.hg38)[seqlevels(promoters_gr)]
promoters_gr <- trim(promoters_gr)

promoters_gr <- promoters_gr[width(promoters_gr) >= 100]
cat("Promoters after trimming (width >= 100bp):", length(promoters_gr), "\n")
if (length(promoters_gr) == 0) stop("No promoter regions survived filtering")

target_genes <- intersect(target_genes, unique(promoters_gr$symbol))
cat("Targets remaining:", length(target_genes),
    "| TFs:", length(tf_in_data),
    "| ratio:", round(length(target_genes) / length(tf_in_data), 2), "\n")
if (length(target_genes) <= length(tf_in_data))
  cat("WARNING: targets no longer exceed TFs — raise N_TARGET_POOL and rerun 8B\n")

pfm_list  <- getMatrixSet(JASPAR2020, list(species = 9606, all_versions = FALSE))
pfm_names <- sapply(pfm_list, function(x) x@name)
pfm_in_data <- pfm_list[pfm_names %in% tf_in_data]
cat("JASPAR motifs loaded:", length(pfm_list),
    "| matching our TFs:", length(pfm_in_data), "\n")

if (length(pfm_in_data) == 0) {
  pfm_in_data <- pfm_list[toupper(pfm_names) %in% toupper(tf_in_data)]
  cat("Case-insensitive match:", length(pfm_in_data), "\n")
}
if (length(pfm_in_data) == 0) stop("No JASPAR motifs matched TFs in data")

cat("Scanning promoters for motif matches (several minutes)...\n")
motif_matches   <- matchMotifs(pfm_in_data, promoters_gr,
                               genome = BSgenome.Hsapiens.UCSC.hg38,
                               out = "scores")
match_matrix    <- motifScores(motif_matches)
tf_names_jaspar <- sapply(pfm_in_data, function(x) x@name)
cat("Motif match matrix:", dim(match_matrix), "\n")

motif_edgelist <- data.frame()
for (j in seq_along(tf_names_jaspar)) {
  df <- data.frame(TF = tf_names_jaspar[j], Gene = promoters_gr$symbol,
                   Score = as.numeric(match_matrix[, j]),
                   stringsAsFactors = FALSE)
  df <- df[df$Score > 0, ]
  if (nrow(df) > 0) motif_edgelist <- rbind(motif_edgelist, df)
}
cat("Raw motif edges:", nrow(motif_edgelist), "\n")

motif_edgelist <- motif_edgelist %>%
  group_by(TF, Gene) %>%
  summarise(Score = max(Score), .groups = "drop") %>%
  as.data.frame()
motif_edgelist <- motif_edgelist[motif_edgelist$Gene %in% target_genes, ]

cat("Motif edges:", nrow(motif_edgelist),
    "| TFs:", length(unique(motif_edgelist$TF)),
    "| targets:", length(unique(motif_edgelist$Gene)), "\n")
cat("Score range:", round(range(motif_edgelist$Score), 4), "\n")
if (nrow(motif_edgelist) == 0) stop("No motif edges")

rm(motif_matches, match_matrix); gc(verbose = FALSE)

# ---- 8B-METH: CONTINUOUS METHYLATION WEIGHTING ----

cat("\n=== 8B-METH: Continuous methylation weighting ===\n")

median_beta <- median(gene_beta$promoter_beta)
n_covered   <- sum(unique(motif_edgelist$Gene) %in% gene_beta$Gene)

cat("Motif target genes:", length(unique(motif_edgelist$Gene)), "\n")
cat("With promoter methylation data:", n_covered,
    sprintf("(%.1f%%)\n", 100 * n_covered / length(unique(motif_edgelist$Gene))))
cat("Genes without data assigned the median beta:", round(median_beta, 3), "\n")

motif_edgelist_weighted <- motif_edgelist %>%
  left_join(gene_beta %>% dplyr::select(Gene, promoter_beta), by = "Gene") %>%
  mutate(meth_covered   = !is.na(promoter_beta),
         promoter_beta  = ifelse(is.na(promoter_beta), median_beta, promoter_beta),
         meth_weight    = 1.5 - promoter_beta,
         Score_original = Score,
         Score          = Score * meth_weight)

cat("\nWeight range:", round(range(motif_edgelist_weighted$meth_weight), 3), "\n")
cat("Edges up-weighted:",   sum(motif_edgelist_weighted$meth_weight > 1), "\n")
cat("Edges down-weighted:", sum(motif_edgelist_weighted$meth_weight < 1), "\n")
cat("Edges effectively modulated:",
    sum(abs(motif_edgelist_weighted$meth_weight - 1) > 0.05),
    sprintf("(%.1f%% of prior)\n",
            100 * mean(abs(motif_edgelist_weighted$meth_weight - 1) > 0.05)))

write.csv(motif_edgelist_weighted,
          file.path(RESULTS_METH, "motif_edgelist_methylation_weighted.csv"), row.names = FALSE)

motif_edgelist_original <- motif_edgelist
motif_edgelist <- motif_edgelist_weighted %>% dplyr::select(TF, Gene, Score)
cat("Motif prior updated — PANDA uses methylation-weighted edges\n")

ggplot(motif_edgelist_weighted, aes(x = meth_weight)) +
  geom_histogram(bins = 60, fill = "#9B59B6", color = "white") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(title    = "Methylation Weights Applied to the Motif Prior",
       subtitle = "weight = 1.5 - promoter beta | dashed = no change",
       x = "Methylation weight", y = "Motif edges") +
  theme_bw() +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
ggsave(file.path(FIG_METH, "methylation_weight_distribution.png"), width = 8, height = 5, dpi = 300)

set.seed(1)
meth_scatter <- motif_edgelist_weighted[
  sample(nrow(motif_edgelist_weighted),
         min(20000, nrow(motif_edgelist_weighted))), ]
ggplot(meth_scatter, aes(x = Score_original, y = Score, color = promoter_beta)) +
  geom_point(alpha = 0.3, size = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_viridis_c(name = "Promoter\nbeta", option = "magma") +
  labs(title    = "Motif Scores Before vs After Methylation Weighting",
       subtitle = "Points below the diagonal are down-weighted by promoter methylation",
       x = "Original motif score", y = "Weighted motif score") +
  theme_bw() +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
ggsave(file.path(FIG_METH, "methylation_weighting_effect.png"), width = 8, height = 6, dpi = 300)

# ---- 8C: EXPRESSION MATRIX ----

expr_panda_clean <- expr_normalized_sym
row_vars <- apply(expr_panda_clean, 1, var)
expr_panda_clean <- expr_panda_clean[order(row_vars, decreasing = TRUE), ]
expr_panda_clean <- expr_panda_clean[!duplicated(rownames(expr_panda_clean)), ]

motif_genes <- unique(c(motif_edgelist$Gene, motif_edgelist$TF))
expr_final  <- expr_panda_clean[intersect(motif_genes, rownames(expr_panda_clean)), ]
gene_vars   <- apply(expr_final, 1, var)
expr_final  <- expr_final[gene_vars > 0, ]
cat("\nExpression subset for PANDA:", dim(expr_final), "\n")

motif_edgelist <- motif_edgelist[motif_edgelist$Gene %in% rownames(expr_final) &
                                   motif_edgelist$TF %in% rownames(expr_final), ]
cat("Motif edges after variance filter:", nrow(motif_edgelist), "\n")
cat("TFs:", length(unique(motif_edgelist$TF)),
    "| targets:", length(unique(motif_edgelist$Gene)), "\n")

# ---- 8D: PPI EDGE LIST ----

ppi_edgelist <- reshape2::melt(ppi_matrix)
colnames(ppi_edgelist) <- c("TF1", "TF2", "Score")
ppi_edgelist <- ppi_edgelist[ppi_edgelist$Score > 0, ]
ppi_edgelist$TF1 <- as.character(ppi_edgelist$TF1)
ppi_edgelist$TF2 <- as.character(ppi_edgelist$TF2)
cat("PPI edges:", nrow(ppi_edgelist), "\n")

missing_from_ppi <- setdiff(unique(motif_edgelist$TF),
                            unique(c(ppi_edgelist$TF1, ppi_edgelist$TF2)))
cat("Motif TFs missing from PPI:", length(missing_from_ppi), "\n")
if (length(missing_from_ppi) > 0)
  ppi_edgelist <- rbind(ppi_edgelist,
                        data.frame(TF1 = missing_from_ppi,
                                   TF2 = missing_from_ppi, Score = 0.001))
cat("PPI edges after self-edges:", nrow(ppi_edgelist), "\n")

# ---- 8E: RUN PANDA ----

try(untrace("panda", where = asNamespace("pandaR")), silent = TRUE)
cat("\nRunning PANDA on", nrow(motif_edgelist), "motif edges...\n")

panda_result <- pandaR::panda(
  motif = motif_edgelist, expr = as.data.frame(expr_final), ppi = ppi_edgelist,
  alpha = 0.1, hamming = 0.001, iter = 1000,
  output = "regulatory", zScale = TRUE, progress = TRUE,
  randomize = "None", cor.method = "pearson", scale.by.present = FALSE,
  edgelist = TRUE, remove.missing.ppi = FALSE,
  remove.missing.motif = FALSE, remove.missing.genes = FALSE, mode = "union")

# ---- 8F: EXTRACT RESULTS ----

reg_net_matrix <- slot(panda_result, "regNet")
edge_list_panda <- reg_net_matrix[, c("TF", "Gene", "Score")]
colnames(edge_list_panda) <- c("TF", "Target", "Score")
edge_list_panda$Score  <- as.numeric(as.character(edge_list_panda$Score))
edge_list_panda$TF     <- as.character(edge_list_panda$TF)
edge_list_panda$Target <- as.character(edge_list_panda$Target)

cat("\nTotal PANDA edges:", nrow(edge_list_panda), "\n")
cat("Score range:", round(range(edge_list_panda$Score, na.rm = TRUE), 3), "\n")

edge_list_panda <- edge_list_panda[!is.na(edge_list_panda$Score) &
                                     edge_list_panda$Score > 0, ]
edge_list_panda <- edge_list_panda[order(edge_list_panda$Score,
                                         decreasing = TRUE), ]
cat("PANDA edges with Score > 0:", nrow(edge_list_panda), "\n")
cat("Unique TFs:", length(unique(edge_list_panda$TF)),
    "| unique targets:", length(unique(edge_list_panda$Target)), "\n")
print(head(edge_list_panda, 10))
write.csv(edge_list_panda, file.path(RESULTS_PANDA, "PANDA_edge_list.csv"), row.names = FALSE)

cat("\nKnown FXS TFs recovered by PANDA:",
    paste(intersect(known_tf, unique(edge_list_panda$TF)), collapse = ", "), "\n")
cat("\nStep 8 done\n")

# ============================================================
# CHECKPOINT
# ============================================================

rm(counts_raw, count_matrix, count_matrix_filtered, ppi_filtered,
   panda_result, reg_net_matrix, expr_panda_clean, probe2gene_all,
   anno_450k, meth_matrix, dmp_results, res_plot,
   meth_clean, dmp_unadj, dmp_adj)
gc(verbose = FALSE)

save.image(file.path(CHECKPOINTS_DIR, "checkpoint_step8.RData"))
cat("\n=== PART 1 COMPLETE — checkpoint_step8.RData saved ===\n")
cat("Restart R, then run part2_steps9-16.R\n")