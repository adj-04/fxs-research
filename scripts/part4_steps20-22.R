# ============================================================
# FXS PIPELINE v5 — PART 4 : STEPS 20–22 (new in v5)
#   20  Tissue-generality classification (paper Sec. IV-J / V-F;
#       in v4 this was computed outside the scripts)
#   21  CollecTRI: curated-regulon validation of inferred targets and
#       TF-activity inference in both cohorts (TF-level replication)
#   22  One summary table of the headline numbers for the manuscript
#
# Requires: checkpoint_step19.RData (from part 3)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggrepel)
  library(org.Hs.eg.db); library(AnnotationDbi)
})

source(file.path(Sys.getenv("FXS_REPO", "."), "scripts", "00_config.R"))
load(ckpt("checkpoint_step19.RData"))
source(file.path(Sys.getenv("FXS_REPO", "."), "scripts", "00_config.R"))
log_session("part4")
cat("Checkpoint loaded.\n")

# ============================================================
# STEP 20 — TISSUE-GENERALITY CLASSIFICATION
# ============================================================

cat("\n=== STEP 20: Tissue-general vs tissue-specific ===\n")

# Current HGNC symbols for every set (450K annotation is hg19-era, e.g.
# ASFMR1 is now FMR1-AS1)
to_current <- function(sym) {
  sym <- unique(sym[!is.na(sym) & sym != ""])
  ok  <- sym %in% keys(org.Hs.eg.db, keytype = "SYMBOL")
  out <- sym
  if (any(!ok)) {
    al <- suppressMessages(mapIds(org.Hs.eg.db, keys = sym[!ok], column = "SYMBOL",
                                  keytype = "ALIAS", multiVals = "first"))
    out[!ok] <- ifelse(is.na(al), sym[!ok], al)
  }
  unique(out)
}

set_org  <- to_current(sig_degs$Symbol)
set_npc  <- to_current(sig_valid$Symbol)
set_meth <- to_current(gene_meth_status$Gene)
# Paper rule: probes in the shared FMR1/ASFMR1 promoter count toward FMR1
if ("FMR1-AS1" %in% set_meth) set_meth <- union(set_meth, "FMR1")

classify <- function(sets) {
  all_genes <- Reduce(union, sets)
  m <- sapply(sets, function(s) all_genes %in% s)
  rownames(m) <- all_genes
  df <- data.frame(Gene = all_genes, m, N_contexts = rowSums(m), check.names = FALSE)
  df$Class <- ifelse(df$N_contexts >= 2, "Tissue-general", "Tissue-specific")
  df[order(-df$N_contexts, df$Gene), ]
}

tissue_df <- classify(list(Organoid = set_org, NPC_B_iPSC = set_npc, Blood_meth = set_meth))
write.csv(tissue_df, file.path(RESULTS_TISSUE, "tissue_generality_classification.csv"),
          row.names = FALSE)

tissue_summary <- data.frame(
  Category = c("Tissue-general (>=2 of 3)", "All three contexts",
               "Specific: organoid only", "Specific: NPC/B/iPSC only",
               "Specific: blood methylation only"),
  N = c(sum(tissue_df$N_contexts >= 2), sum(tissue_df$N_contexts == 3),
        sum(tissue_df$N_contexts == 1 & tissue_df$Organoid),
        sum(tissue_df$N_contexts == 1 & tissue_df$NPC_B_iPSC),
        sum(tissue_df$N_contexts == 1 & tissue_df$Blood_meth)),
  Genes = c(paste(tissue_df$Gene[tissue_df$N_contexts >= 2], collapse = ";"),
            paste(tissue_df$Gene[tissue_df$N_contexts == 3], collapse = ";"),
            "", "", ""))
print(tissue_summary[, 1:2])
write.csv(tissue_summary, file.path(RESULTS_TISSUE, "tissue_generality_summary.csv"),
          row.names = FALSE)

# Chance expectation for each pairwise overlap (hypergeometric on the
# shared background) so "only FMR1" can be read against what chance gives
bg <- to_current(intersect(res_df$Symbol, deg_valid$Symbol))
pair_tests <- do.call(rbind, combn(c("Organoid", "NPC_B_iPSC", "Blood_meth"), 2,
  simplify = FALSE, FUN = function(pr) {
    sets <- list(Organoid = set_org, NPC_B_iPSC = set_npc, Blood_meth = set_meth)
    a <- intersect(sets[[pr[1]]], bg); b <- intersect(sets[[pr[2]]], bg)
    k <- length(intersect(a, b)); N <- length(bg)
    data.frame(Pair = paste(pr, collapse = " vs "), Size_A = length(a), Size_B = length(b),
               Overlap = k, Expected = length(a) * length(b) / N,
               p_hyper = phyper(k - 1, length(b), N - length(b), length(a), lower.tail = FALSE),
               Genes = paste(intersect(a, b), collapse = ";"))
  }))
print(pair_tests)
write.csv(pair_tests, file.path(RESULTS_TISSUE, "tissue_pairwise_overlap_tests.csv"), row.names = FALSE)

# Sensitivity: validation DEGs from the donor-blocked limma model
if (exists("sig_valid_limma")) {
  tissue_df_l <- classify(list(Organoid = set_org,
                               NPC_B_iPSC = to_current(sig_valid_limma$Symbol),
                               Blood_meth = set_meth))
  write.csv(tissue_df_l, file.path(RESULTS_TISSUE,
            "tissue_generality_classification_donorblocked.csv"), row.names = FALSE)
  cat("Tissue-general genes with donor-blocked validation DEGs:",
      paste(tissue_df_l$Gene[tissue_df_l$N_contexts >= 2], collapse = ", "), "\n")
}

ggplot(tissue_summary[c(1, 3:5), ], aes(x = reorder(Category, N), y = N,
                                         fill = grepl("general", Category))) +
  geom_col() + coord_flip() +
  geom_text(aes(label = N), hjust = -0.2) +
  scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#95A5A6"), guide = "none") +
  labs(title = "Tissue-general vs tissue-specific significance",
       subtitle = "Organoid RNA-seq | NPC/B/iPSC RNA-seq | blood 450K methylation",
       x = NULL, y = "Genes") +
  theme_bw() + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(FIG_TISSUE, "tissue_generality.png"), width = 8, height = 4.5, dpi = 300)

# ============================================================
# STEP 21 — CollecTRI VALIDATION + TF ACTIVITY
# ============================================================

cat("\n=== STEP 21: CollecTRI ===\n")
collectri <- tryCatch({
  suppressPackageStartupMessages(library(decoupleR))
  decoupleR::get_collectri(organism = "human", split_complexes = FALSE)
}, error = function(e) { cat("CollecTRI unavailable:", conditionMessage(e), "\n"); NULL })

if (!is.null(collectri)) {
  collectri <- as.data.frame(collectri)
  write.csv(data.frame(Interactions = nrow(collectri),
                       TFs = length(unique(collectri$source)),
                       Retrieved = format(Sys.Date())),
            file.path(RESULTS_VALID, "CollecTRI_version.csv"), row.names = FALSE)

  # ---- 21A: curated-regulon enrichment of predicted targets ----
  universe_genes <- unique(c(edge_list_genie3$Target, edge_list_panda$Target,
                             edge_list_genie3$TF, edge_list_panda$TF))
  test_tfs <- intersect(validation_tfs, unique(collectri$source))
  collectri_val <- do.call(rbind, lapply(test_tfs, function(tf) {
    cur <- intersect(unique(collectri$target[collectri$source == tf]), universe_genes)
    if (length(cur) < 3) return(NULL)
    do.call(rbind, lapply(c("GENIE3", "PANDA"), function(src) {
      pred <- intersect(unique(get_tf_targets(tf, src, top_n = 500)), universe_genes)
      if (length(pred) < 3) return(NULL)
      ov <- length(intersect(pred, cur)); ex <- length(pred) * length(cur) / length(universe_genes)
      data.frame(TF = tf, Source = src, N_predicted = length(pred), N_curated = length(cur),
                 N_overlap = ov, Expected = ex, Fold_enrich = ov / ex,
                 pvalue = phyper(ov - 1, length(cur), length(universe_genes) - length(cur),
                                 length(pred), lower.tail = FALSE))
    }))
  }))
  if (!is.null(collectri_val) && nrow(collectri_val) > 0) {
    collectri_val$FDR <- p.adjust(collectri_val$pvalue, method = "BH")
    collectri_val <- collectri_val[order(collectri_val$pvalue), ]
    print(head(collectri_val, 20))
    write.csv(collectri_val, file.path(RESULTS_VALID, "CollecTRI_validation_results.csv"),
              row.names = FALSE)
  }

  # ---- 21B: TF activity from each cohort's DE statistics ----
  # Activity is inferred from a TF's whole curated regulon, so it can
  # replicate across cohorts even when individual target genes do not.
  stat_mat <- function(df) {
    d <- df %>% filter(!is.na(stat), !is.na(Symbol), Symbol != "") %>%
      group_by(Symbol) %>% slice_max(abs(stat), n = 1, with_ties = FALSE) %>% ungroup()
    matrix(d$stat, ncol = 1, dimnames = list(d$Symbol, "stat"))
  }
  act <- tryCatch(lapply(list(Organoid = res_df, NPC_B_iPSC = deg_valid), function(df) {
    r <- decoupleR::run_ulm(mat = stat_mat(df), network = collectri, .source = "source",
                            .target = "target", .mor = "mor", minsize = 5)
    as.data.frame(r)[, c("source", "score", "p_value")]
  }), error = function(e) { cat("TF activity failed:", conditionMessage(e), "\n"); NULL })
}
if (!is.null(collectri) && !is.null(act)) {
  tf_act <- merge(act$Organoid, act$NPC_B_iPSC, by = "source",
                  suffixes = c("_organoid", "_npc"))
  tf_act$FDR_organoid <- p.adjust(tf_act$p_value_organoid, "BH")
  tf_act$FDR_npc      <- p.adjust(tf_act$p_value_npc, "BH")
  tf_act$Known_FXS    <- tf_act$source %in% known_tf
  tf_act <- tf_act[order(tf_act$p_value_organoid), ]
  write.csv(tf_act, file.path(RESULTS_VALID, "TF_activity_CollecTRI_both_cohorts.csv"),
            row.names = FALSE)

  act_cor <- cor.test(tf_act$score_organoid, tf_act$score_npc, method = "spearman", exact = FALSE)
  both_sig <- tf_act %>% filter(FDR_organoid < 0.05, FDR_npc < 0.05)
  act_conc <- mean(sign(both_sig$score_organoid) == sign(both_sig$score_npc))
  cat("TF activity cross-cohort Spearman rho:", round(act_cor$estimate, 3),
      "p =", signif(act_cor$p.value, 3), "| TFs:", nrow(tf_act), "\n")
  cat("TFs FDR<0.05 in both cohorts:", nrow(both_sig),
      "| same direction:", round(100 * act_conc, 1), "%\n")
  cat("Known FXS TFs:\n")
  print(tf_act[tf_act$Known_FXS, c("source", "score_organoid", "FDR_organoid",
                                   "score_npc", "FDR_npc")])
  write.csv(data.frame(Metric = c("TFs scored in both", "Spearman rho", "Spearman p",
                                  "TFs FDR<0.05 in both", "Same-direction (%)"),
                       Value = c(nrow(tf_act), round(act_cor$estimate, 4), act_cor$p.value,
                                 nrow(both_sig), round(100 * act_conc, 1))),
            file.path(RESULTS_VALID, "TF_activity_crosscohort_summary.csv"), row.names = FALSE)

  ggplot(tf_act, aes(score_organoid, score_npc)) +
    geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
    geom_point(alpha = 0.4, colour = "grey50") +
    geom_point(data = subset(tf_act, Known_FXS), colour = "#E74C3C", size = 2.5) +
    geom_text_repel(data = subset(tf_act, Known_FXS | (FDR_organoid < 0.05 & FDR_npc < 0.05)),
                    aes(label = source), size = 3, max.overlaps = 25) +
    labs(title = "TF activity (CollecTRI, ULM) — organoid vs NPC/B/iPSC",
         subtitle = sprintf("Spearman rho = %.3f, p = %.2g | red = a priori FXS TFs",
                            act_cor$estimate, act_cor$p.value),
         x = "Activity score, GSE146339", y = "Activity score, GSE218676") +
    theme_bw() + theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(FIG_VALID, "TF_activity_crosscohort.png"), width = 8, height = 7, dpi = 300)
}

# ============================================================
# STEP 22 — HEADLINE NUMBERS FOR THE MANUSCRIPT
# ============================================================

cat("\n=== STEP 22: Summary for the manuscript ===\n")
rd <- function(f) if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL
getv <- function(df, m) if (is.null(df)) NA else df$Value[match(m, df$Metric)]

g3n <- rd(file.path(RESULTS_GENIE3, "GENIE3_null_summary.csv"))
pdn <- rd(file.path(RESULTS_PANDA, "PANDA_null_summary.csv"))
cc  <- rd(file.path(RESULTS_VALID, "crossvalidation_concordance_summary.csv"))
ds  <- rd(file.path(RESULTS_VALID, "validation_donor_sensitivity.csv"))
mm  <- rd(file.path(RESULTS_METH, "DMP_model_comparison.csv"))
asg <- rd(file.path(RESULTS_ALPACA, "ALPACA_results_with_significance.csv"))
msg <- rd(file.path(RESULTS_ALPACA, "ALPACA_module_significance.csv"))

summary_tbl <- data.frame(
  Item = c("Primary DEGs", "Validation DEGs (DESeq2)", "Validation DEGs (limma, donor-blocked)",
           "Sig. DMPs (cell-adjusted)", "Sig. DMPs (unadjusted, v4 model)",
           "GENIE3 genes / non-TF targets", "GENIE3 edges pooled FDR<0.05",
           "PANDA edges pooled FDR<0.05", "ALPACA nodes pooled FDR<0.05",
           "Best ALPACA module FDR", "LIONESS edges limma FDR<0.05 (timepoint-adjusted)",
           "Cross-cohort Pearson r (raw LFC)", "Cross-cohort Pearson r (shrunken LFC)",
           "RRHO peak -log10 p", "Tissue-general genes"),
  Value = c(nrow(sig_degs), nrow(sig_valid), getv(ds, "DEGs limma-voom (donor blocked)"),
            if (!is.null(mm)) mm$Sig_DMPs[2] else NA, if (!is.null(mm)) mm$Sig_DMPs[1] else NA,
            paste(nrow(expr_subset), length(setdiff(rownames(expr_subset), tf_in_data)), sep = " / "),
            getv(g3n, "Edges pooled FDR < 0.05"), getv(pdn, "Edges pooled FDR < 0.05"),
            if (!is.null(asg)) sum(asg$emp_FDR < 0.05, na.rm = TRUE) else NA,
            if (!is.null(msg)) signif(min(msg$module_FDR), 3) else NA,
            sum(diff_df$limma_FDR < 0.05),
            getv(cc, "Pearson r"), getv(cc, "Pearson r (shrunken LFC)"),
            getv(cc, "RRHO peak -log10(p)"),
            paste(tissue_df$Gene[tissue_df$N_contexts >= 2], collapse = ", ")))
# v5.1: the v5 ALPACA node/module tests are invalid (module scores sum to 2);
# drop them and add the robust-null rows from part 5 when available.
summary_tbl$Value <- as.character(summary_tbl$Value)
summary_tbl <- summary_tbl[!summary_tbl$Item %in% c("ALPACA nodes pooled FDR<0.05", "Best ALPACA module FDR"), ]
summary_tbl <- rbind(summary_tbl, robust_headline_rows())
print(summary_tbl, row.names = FALSE)
write.csv(summary_tbl, file.path("results", "v5_headline_numbers.csv"), row.names = FALSE)

save.image(ckpt("checkpoint_step22.RData"))
cat("\n=== PIPELINE v5 COMPLETE ===\n")
