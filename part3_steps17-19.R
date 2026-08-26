# ============================================================
# FXS PIPELINE v4 — PART 3 of 3 : STEPS 17–19
# Permutation nulls, weighted centrality, in silico validation
#
# Requires: checkpoint_step16.RData (from part 2)
#
# MEMORY NOTES (this is the part that crashes on low-RAM machines)
#   - N_CORES = 2: GENIE3 forks the expression matrix per core
#   - gc() inside the permutation loop: R otherwise holds discarded
#     matrices until the next allocation forces collection
#   - Step 17A attaches empirical p by matrix indexing rather than
#     melting a multi-million-row matrix and left_join-ing it, which
#     transiently held three copies
#   - RUN_PANDA_NULL defaults FALSE: run it separately once the rest
#     of Step 17 completes
# ============================================================

library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(GENIE3)
library(pandaR)
library(netZooR)
library(igraph)
library(httr)
library(jsonlite)

setwd("~/Desktop/final report")
load("checkpoint_step16.RData")
cat("Checkpoint loaded. Objects:", length(ls()), "\n")

# ------------------------------------------------------------
# CONFIG — tune down first, raise once it completes
# ------------------------------------------------------------

N_PERM_GENIE3  <- 50      # 100 if memory allows
N_TREES_PERM   <- 250
N_PERM_PANDA   <- 20
N_PERM_ALPACA  <- 50
N_BOOT_LIONESS <- 1000
N_CORES        <- 2       # forking multiplies memory; keep low
RUN_PANDA_NULL <- FALSE   # TRUE only after 17A/C/D succeed
set.seed(42)

cat("\n=== STEP 17: Network significance and robustness ===\n")

# ---- 17A: GENIE3 PERMUTATION NULL ----
# GENIE3 never sees the FXS/Control labels, so permuting labels is a
# no-op. The correct null shuffles each gene's expression vector
# independently, destroying co-expression while preserving marginals.

cat("\n--- 17A: GENIE3 permutation null (", N_PERM_GENIE3, "permutations) ---\n")

permute_rows <- function(m) {
  out <- t(apply(m, 1, sample)); dimnames(out) <- dimnames(m); out
}

exceed_count <- matrix(0L, nrow = nrow(genie3_weights),
                       ncol = ncol(genie3_weights),
                       dimnames = dimnames(genie3_weights))
null_pool <- numeric(0)

t0 <- Sys.time()
for (p in seq_len(N_PERM_GENIE3)) {
  if (p %% 10 == 0 || p == 1)
    cat("  permutation", p, "/", N_PERM_GENIE3, "| elapsed",
        round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
  w_perm <- GENIE3(permute_rows(expr_subset), regulators = tf_in_data,
                   nTrees = N_TREES_PERM, nCores = N_CORES, verbose = FALSE)
  w_perm <- w_perm[rownames(genie3_weights), colnames(genie3_weights),
                   drop = FALSE]
  exceed_count <- exceed_count + (w_perm >= genie3_weights)
  null_pool <- c(null_pool, sample(as.numeric(w_perm),
                                   size = min(50000, length(w_perm))))
  rm(w_perm); gc(verbose = FALSE)          # release before next allocation
}
cat("Complete in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

genie3_emp_p <- (exceed_count + 1) / (N_PERM_GENIE3 + 1)
rm(exceed_count); gc(verbose = FALSE)

null_q95 <- quantile(null_pool, 0.95)
null_q99 <- quantile(null_pool, 0.99)
cat("Null 95th pct:", signif(null_q95, 4),
    "| 99th pct:", signif(null_q99, 4), "\n")
cat("Observed edges above q95:", sum(genie3_weights >= null_q95),
    "| above q99:", sum(genie3_weights >= null_q99), "\n")

# Matrix indexing instead of melt + left_join. The old version built a
# data frame with one row per matrix cell and held two copies during
# the join — the single largest allocation in the pipeline.
idx <- cbind(match(as.character(edge_list_genie3$Target), rownames(genie3_emp_p)),
             match(as.character(edge_list_genie3$TF),     colnames(genie3_emp_p)))

edge_list_genie3_sig <- edge_list_genie3
edge_list_genie3_sig$emp_p        <- genie3_emp_p[idx]
edge_list_genie3_sig$emp_FDR      <- p.adjust(edge_list_genie3_sig$emp_p,
                                              method = "BH")
edge_list_genie3_sig$above_null95 <- edge_list_genie3_sig$Weight >= null_q95
edge_list_genie3_sig$above_null99 <- edge_list_genie3_sig$Weight >= null_q99
edge_list_genie3_sig <- edge_list_genie3_sig[
  order(edge_list_genie3_sig$emp_p, -edge_list_genie3_sig$Weight), ]

rm(idx); gc(verbose = FALSE)

cat("Edges with empirical p < 0.05:",
    sum(edge_list_genie3_sig$emp_p < 0.05, na.rm = TRUE), "\n")
cat("Edges with empirical FDR < 0.05:",
    sum(edge_list_genie3_sig$emp_FDR < 0.05, na.rm = TRUE), "\n")
cat("NOTE: with this many edges tested, BH is severe. Report the\n")
cat("nominal and FDR counts honestly rather than only the former.\n")
write.csv(edge_list_genie3_sig, "GENIE3_edge_list_with_null.csv", row.names = FALSE)

null_plot_df <- rbind(
  data.frame(Weight = sample(null_pool, min(100000, length(null_pool))),
             Source = "Permuted null"),
  data.frame(Weight = sample(as.numeric(genie3_weights),
                             min(100000, length(genie3_weights))),
             Source = "Observed"))
null_plot_df <- null_plot_df[null_plot_df$Weight > 0, ]

ggplot(null_plot_df, aes(x = Weight, fill = Source)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = null_q95, linetype = "dashed", color = "#2C3E50") +
  scale_x_log10() +
  scale_fill_manual(values = c("Observed" = "#E74C3C",
                               "Permuted null" = "#95A5A6")) +
  labs(title = "GENIE3 Edge Weights — Observed vs Permuted Null",
       subtitle = paste0(N_PERM_GENIE3,
                         " per-gene permutations | dashed = null 95th pct"),
       x = "Edge weight (log10)", y = "Density") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
ggsave("GENIE3_null_distribution.png", width = 9, height = 6, dpi = 300)

rm(null_plot_df, null_pool); gc(verbose = FALSE)

genie3_null_filtered <- edge_list_genie3_sig %>%
  filter(above_null95, emp_p < 0.05)
cat("GENIE3 edges surviving null filter:", nrow(genie3_null_filtered), "\n")
if (nrow(genie3_null_filtered) > 0) {
  g_genie3_null <- graph_from_data_frame(
    data.frame(from = genie3_null_filtered$TF, to = genie3_null_filtered$Target,
               weight = genie3_null_filtered$Weight), directed = TRUE)
  V(g_genie3_null)$type <- ifelse(V(g_genie3_null)$name %in% tf_in_data,
                                  "TF", "Target")
  write_graph(g_genie3_null, "GENIE3_network_nullfiltered.graphml",
              format = "graphml")
  rm(g_genie3_null)
}
gc(verbose = FALSE)

# ---- 17B: PANDA RANDOMIZED PRIOR ----

if (RUN_PANDA_NULL) {
  cat("\n--- 17B: PANDA randomized-prior null (", N_PERM_PANDA, "runs) ---\n")

  randomize_motif_degpres <- function(motif_df) {
    bg <- graph_from_data_frame(
      data.frame(from = paste0("TF_", motif_df$TF),
                 to   = paste0("G_",  motif_df$Gene)), directed = FALSE)
    bg_rw <- rewire(bg, keeping_degseq(niter = ecount(bg) * 10))
    el <- as_edgelist(bg_rw)
    flip <- grepl("^G_", el[, 1])
    el[flip, ] <- el[flip, c(2, 1)]
    data.frame(TF = sub("^TF_", "", el[, 1]), Gene = sub("^G_", "", el[, 2]),
               Score = motif_df$Score[sample(nrow(motif_df))],
               stringsAsFactors = FALSE)
  }

  obs_scores <- setNames(edge_list_panda$Score,
                         paste(edge_list_panda$TF, edge_list_panda$Target,
                               sep = "__"))
  panda_exceed <- setNames(rep(0L, length(obs_scores)), names(obs_scores))
  panda_null_pool <- numeric(0)

  for (p in seq_len(N_PERM_PANDA)) {
    cat("  PANDA null run", p, "/", N_PERM_PANDA, "\n")
    motif_rand <- randomize_motif_degpres(motif_edgelist)
    motif_rand <- motif_rand[motif_rand$Gene %in% rownames(expr_final) &
                               motif_rand$TF %in% rownames(expr_final), ]
    if (nrow(motif_rand) == 0) next
    pr <- tryCatch(
      pandaR::panda(motif = motif_rand, expr = as.data.frame(expr_final),
                    ppi = ppi_edgelist, alpha = 0.1, hamming = 0.001,
                    iter = 1000, output = "regulatory", zScale = TRUE,
                    progress = FALSE, randomize = "None",
                    cor.method = "pearson", scale.by.present = FALSE,
                    edgelist = TRUE, remove.missing.ppi = FALSE,
                    remove.missing.motif = FALSE, remove.missing.genes = FALSE,
                    mode = "union"),
      error = function(e) { cat("    failed:", conditionMessage(e), "\n"); NULL })
    if (is.null(pr)) next
    rn <- slot(pr, "regNet")
    null_edges <- setNames(as.numeric(as.character(rn$Score)),
                           paste(as.character(rn$TF), as.character(rn$Gene),
                                 sep = "__"))
    panda_null_pool <- c(panda_null_pool,
                         sample(null_edges, min(50000, length(null_edges))))
    common <- intersect(names(obs_scores), names(null_edges))
    panda_exceed[common] <- panda_exceed[common] +
      (null_edges[common] >= obs_scores[common])
    rm(pr, rn, null_edges); gc(verbose = FALSE)
  }

  panda_emp_p <- (panda_exceed + 1) / (N_PERM_PANDA + 1)
  panda_null_q95 <- quantile(panda_null_pool, 0.95, na.rm = TRUE)

  edge_list_panda_sig <- edge_list_panda %>%
    mutate(Edge = paste(TF, Target, sep = "__"),
           emp_p = panda_emp_p[Edge],
           emp_FDR = p.adjust(panda_emp_p[Edge], method = "BH"),
           above_null95 = Score >= panda_null_q95) %>%
    arrange(emp_p, desc(Score))

  cat("PANDA null 95th pct:", signif(panda_null_q95, 4), "\n")
  cat("Observed 90th pct was:", signif(panda_threshold, 4), "\n")
  cat("Edges above randomized-prior null q95:",
      sum(edge_list_panda_sig$above_null95, na.rm = TRUE), "\n")
  cat("INTERPRET: if the null q95 sits close to the observed 90th\n")
  cat("percentile, PANDA is largely reproducing its prior rather than\n")
  cat("adding expression-driven signal. Report this in Limitations.\n")
  write.csv(edge_list_panda_sig, "PANDA_edge_list_with_null.csv", row.names = FALSE)

  panda_cmp <- rbind(
    data.frame(Score = sample(as.numeric(panda_null_pool),
                              min(100000, length(panda_null_pool))),
               Source = "Randomized prior"),
    data.frame(Score = sample(edge_list_panda$Score,
                              min(100000, nrow(edge_list_panda))),
               Source = "Observed"))
  ggplot(panda_cmp, aes(x = Score, fill = Source)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = panda_null_q95, linetype = "dashed") +
    scale_fill_manual(values = c("Observed" = "#E74C3C",
                                 "Randomized prior" = "#95A5A6")) +
    labs(title = "PANDA Edge Scores — Observed vs Randomized Prior",
         subtitle = paste0(N_PERM_PANDA, " degree-preserving randomizations"),
         x = "PANDA edge score (z-scaled)", y = "Density") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  ggsave("PANDA_null_distribution.png", width = 9, height = 6, dpi = 300)
  rm(panda_cmp, panda_null_pool); gc(verbose = FALSE)
} else {
  cat("\n--- 17B skipped (RUN_PANDA_NULL = FALSE) ---\n")
  cat("Run it separately later: set TRUE, reload the checkpoint, run 17B only.\n")
}

# ---- 17C: ALPACA MODULE SIGNIFICANCE ----
# CRANE output is not a valid pandaToAlpaca input in this netZooR build,
# so weight permutation over fixed topology is used. A probe run
# determines the actual output filename before the loop.

cat("\n--- 17C: ALPACA module significance (", N_PERM_ALPACA, "runs) ---\n")

cat("Probe run to determine pandaToAlpaca output filename...\n")
probe_stem <- "./tmp_alpaca_probe"
probe_net  <- fxs_net
probe_net$Score <- fxs_net$Score[sample(nrow(fxs_net))]
probe <- tryCatch(
  netZooR::pandaToAlpaca(panda.net1 = probe_net, panda.net2 = control_net,
                         file.stem = probe_stem, verbose = FALSE),
  error = function(e) { cat("PROBE FAILED:", conditionMessage(e), "\n"); NULL })

probe_files <- list.files(".", pattern = "tmp_alpaca_probe")
cat("Files written by probe:\n"); print(probe_files)

score_file_name <- grep("scores", probe_files, value = TRUE)
if (length(score_file_name) == 0) {
  cat("WARNING: no scores file from probe. ALPACA null will be skipped.\n")
  score_suffix <- NA
} else {
  score_suffix <- sub("^tmp_alpaca_probe", "", score_file_name[1])
  cat("Score file suffix:", score_suffix, "\n")
}
if (length(probe_files) > 0) file.remove(probe_files)

alpaca_null_scores <- list()
n_failed <- 0

if (!is.na(score_suffix)) {
  for (p in seq_len(N_PERM_ALPACA)) {
    if (p %% 10 == 0 || p == 1)
      cat("  ALPACA null run", p, "/", N_PERM_ALPACA, "\n")

    fxs_null <- fxs_net
    fxs_null$Score <- fxs_net$Score[sample(nrow(fxs_net))]

    stem <- paste0("./tmp_alpaca_null_", p)
    ar <- tryCatch(
      netZooR::pandaToAlpaca(panda.net1 = fxs_null, panda.net2 = control_net,
                             file.stem = stem, verbose = FALSE),
      error = function(e) {
        n_failed <<- n_failed + 1
        if (n_failed <= 3) cat("    run", p, "failed:", conditionMessage(e), "\n")
        NULL })

    sf <- paste0(stem, score_suffix)
    if (file.exists(sf)) {
      sc <- read.table(sf, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                       col.names = c("Node", "Score"))
      alpaca_null_scores[[length(alpaca_null_scores) + 1]] <- sc
    }
    file.remove(list.files(".", pattern = paste0("^tmp_alpaca_null_", p, "[_.]"),
                           full.names = TRUE))
    rm(ar, fxs_null); gc(verbose = FALSE)
  }
}

cat("Successful null runs:", length(alpaca_null_scores), "/", N_PERM_ALPACA, "\n")

if (length(alpaca_null_scores) >= 5) {
  null_all   <- do.call(rbind, alpaca_null_scores)
  null_split <- split(null_all$Score, null_all$Node)

  node_null <- null_all %>%
    group_by(Node) %>%
    summarise(null_mean = mean(Score, na.rm = TRUE),
              null_sd = sd(Score, na.rm = TRUE),
              n_null = n(), .groups = "drop")

  alpaca_sig <- alpaca_df %>% left_join(node_null, by = "Node")
  alpaca_sig$emp_p <- mapply(function(node, score) {
    nd <- null_split[[node]]
    if (is.null(nd) || length(nd) == 0) return(NA_real_)
    (sum(nd >= score) + 1) / (length(nd) + 1)
  }, alpaca_sig$Node, alpaca_sig$Score)
  alpaca_sig$emp_FDR <- p.adjust(alpaca_sig$emp_p, method = "BH")
  alpaca_sig$z_vs_null <- (alpaca_sig$Score - alpaca_sig$null_mean) /
    ifelse(is.na(alpaca_sig$null_sd) | alpaca_sig$null_sd == 0, NA,
           alpaca_sig$null_sd)

  cat("ALPACA nodes with empirical p < 0.05:",
      sum(alpaca_sig$emp_p < 0.05, na.rm = TRUE), "\n")
  cat("ALPACA nodes with empirical FDR < 0.05:",
      sum(alpaca_sig$emp_FDR < 0.05, na.rm = TRUE), "\n")

  module_sig <- alpaca_sig %>%
    group_by(Module) %>%
    summarise(n_nodes = n(), mean_score = mean(Score, na.rm = TRUE),
              mean_null = mean(null_mean, na.rm = TRUE),
              n_sig_nodes = sum(emp_p < 0.05, na.rm = TRUE),
              module_emp_p = {
                obs <- mean(Score, na.rm = TRUE)
                nulls <- replicate(1000, mean(sample(null_all$Score, n()),
                                              na.rm = TRUE))
                (sum(nulls >= obs) + 1) / 1001
              }, .groups = "drop") %>%
    mutate(module_FDR = p.adjust(module_emp_p, method = "BH")) %>%
    arrange(module_emp_p)

  cat("\nModule-level significance:\n"); print(as.data.frame(module_sig))
  write.csv(alpaca_sig, "ALPACA_results_with_significance.csv", row.names = FALSE)
  write.csv(module_sig, "ALPACA_module_significance.csv", row.names = FALSE)

  cat("\nKnown FXS TFs with ALPACA significance:\n")
  print(alpaca_sig[alpaca_sig$Gene %in% known_tf,
                   c("Gene","Module","Score","z_vs_null","emp_p","emp_FDR")])

  ggplot(module_sig %>% filter(n_nodes >= 5),
         aes(x = reorder(factor(Module), -log10(module_emp_p)),
             y = -log10(module_emp_p), fill = module_emp_p < 0.05)) +
    geom_bar(stat = "identity") + coord_flip() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#95A5A6"),
                      name = "p < 0.05") +
    labs(title = "ALPACA Module Significance vs Permuted Null",
         subtitle = paste0(length(alpaca_null_scores),
                           " weight-permuted null networks"),
         x = "Module", y = "-log10(empirical p)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  ggsave("ALPACA_module_significance.png", width = 8, height = 6, dpi = 300)

  rm(null_all, null_split, alpaca_null_scores); gc(verbose = FALSE)
} else {
  cat("Too few successful null runs — do not report ALPACA significance.\n")
}

# ---- 17D: LIONESS BOOTSTRAP ----

cat("\n--- 17D: LIONESS bootstrap (", N_BOOT_LIONESS, "replicates) ---\n")

lioness_mat <- as.matrix(lioness_df)
obs_prop_up <- mean(rowMeans(lioness_mat[, fxs_cols, drop = FALSE]) >
                      rowMeans(lioness_mat[, control_cols, drop = FALSE]))
obs_mean_diff <- mean(rowMeans(lioness_mat[, fxs_cols, drop = FALSE]) -
                        rowMeans(lioness_mat[, control_cols, drop = FALSE]))

boot_stats <- t(replicate(N_BOOT_LIONESS, {
  f  <- sample(fxs_cols,     length(fxs_cols),     replace = TRUE)
  c_ <- sample(control_cols, length(control_cols), replace = TRUE)
  fm <- rowMeans(lioness_mat[, f,  drop = FALSE])
  cm <- rowMeans(lioness_mat[, c_, drop = FALSE])
  c(prop_up = mean(fm > cm), mean_diff = mean(fm - cm))
}))

ci_prop <- quantile(boot_stats[, "prop_up"],   c(0.025, 0.975), na.rm = TRUE)
ci_diff <- quantile(boot_stats[, "mean_diff"], c(0.025, 0.975), na.rm = TRUE)

cat("Edges stronger in FXS:", round(obs_prop_up * 100, 1), "%\n")
cat("  95% CI: [", round(ci_prop[1] * 100, 1), "%, ",
    round(ci_prop[2] * 100, 1), "%]\n", sep = "")
cat("Mean edge difference:", signif(obs_mean_diff, 4), "\n")
cat("  95% CI: [", signif(ci_diff[1], 4), ", ", signif(ci_diff[2], 4), "]\n", sep = "")
if (ci_prop[1] < 0.5 && ci_prop[2] > 0.5)
  cat("CI spans 50% — no detectable directional shift. Report as such.\n")

cat("\nSample sizes: n =", length(fxs_cols), "FXS,",
    length(control_cols), "Control.\n")
cat("Bootstrap CI is indicative at this n — report as a range, not a test.\n")

write.csv(data.frame(
  Statistic = c("Prop. edges stronger in FXS", "Mean edge diff (FXS-Control)"),
  Observed = c(obs_prop_up, obs_mean_diff),
  CI_lower = c(ci_prop[1], ci_diff[1]), CI_upper = c(ci_prop[2], ci_diff[2]),
  N_boot = N_BOOT_LIONESS, N_FXS = length(fxs_cols),
  N_Control = length(control_cols)),
  "LIONESS_bootstrap_CI.csv", row.names = FALSE)

ggplot(data.frame(prop = boot_stats[, "prop_up"]), aes(x = prop)) +
  geom_histogram(bins = 50, fill = "#3498DB", color = "white") +
  geom_vline(xintercept = obs_prop_up, color = "#E74C3C", linewidth = 1) +
  geom_vline(xintercept = ci_prop, linetype = "dashed", color = "#2C3E50") +
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "grey40") +
  labs(title = "LIONESS Bootstrap — Proportion of Edges Stronger in FXS",
       subtitle = paste0("Observed = ", round(obs_prop_up * 100, 1),
                         "% | 95% CI [", round(ci_prop[1] * 100, 1), "%, ",
                         round(ci_prop[2] * 100, 1), "%]"),
       x = "Bootstrap proportion", y = "Count") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
ggsave("LIONESS_bootstrap_distribution.png", width = 8, height = 5, dpi = 300)

rm(lioness_mat, boot_stats); gc(verbose = FALSE)
cat("Step 17 done\n")

# ============================================================
# STEP 18 — WEIGHTED DEGREE + MODULE-LOCAL CENTRALITY
# ============================================================

cat("\n=== STEP 18: Weighted degree and module-local centrality ===\n")

# ---- 18A: WEIGHT-BASED CENTRALITY ----
# GENIE3 normalises each regulator's weights to sum to 1, so summed
# weight carries no information. Mean and max edge weight are used.

cat("\n--- 18A: Raw vs weight-based centrality ---\n")

genie3_strength <- edge_list_genie3 %>%
  group_by(TF) %>%
  summarise(Raw_OutDegree = n(),
            Sum_Weight = sum(Weight),
            Mean_Edge_Weight = mean(Weight),
            Max_Edge_Weight = max(Weight),
            Top1pct_Edges = sum(Weight >= quantile(edge_list_genie3$Weight, 0.99)),
            .groups = "drop") %>%
  mutate(Rank_Raw     = rank(-Raw_OutDegree,    ties.method = "min"),
         Rank_MeanW   = rank(-Mean_Edge_Weight, ties.method = "min"),
         Rank_MaxW    = rank(-Max_Edge_Weight,  ties.method = "min"),
         Rank_Top1pct = rank(-Top1pct_Edges,    ties.method = "min"),
         Rank_Shift   = Rank_Raw - Rank_MeanW) %>%
  arrange(Rank_MeanW)

cat("Sum_Weight range across TFs:",
    round(range(genie3_strength$Sum_Weight), 4), "\n")
cat("(All near 1 — GENIE3 normalises per regulator. Using mean/max.)\n\n")

cat("Top 20 TFs by MEAN edge weight (GENIE3):\n")
print(as.data.frame(head(genie3_strength[, c("TF","Raw_OutDegree",
                                             "Mean_Edge_Weight","Max_Edge_Weight",
                                             "Top1pct_Edges","Rank_Raw",
                                             "Rank_MeanW","Rank_Shift")], 20)))

cat("\nKnown FXS TFs — raw vs weight-based (GENIE3):\n")
print(as.data.frame(genie3_strength %>% filter(TF %in% known_tf) %>%
                      dplyr::select(TF, Raw_OutDegree, Mean_Edge_Weight,
                                    Max_Edge_Weight, Top1pct_Edges,
                                    Rank_Raw, Rank_MeanW, Rank_MaxW, Rank_Shift)))
cat("Total TFs ranked:", nrow(genie3_strength), "\n")
cat("\nFOR DISCUSSION: a TF ranking low on raw degree but high on mean\n")
cat("edge weight regulates few targets with high confidence. Raw degree\n")
cat("systematically hides that pattern.\n")

panda_strength <- edge_list_panda %>%
  group_by(TF) %>%
  summarise(Raw_OutDegree = n(), Weighted_OutDegree = sum(Score),
            Mean_Edge_Weight = mean(Score), Max_Edge_Weight = max(Score),
            .groups = "drop") %>%
  mutate(Rank_Raw = rank(-Raw_OutDegree, ties.method = "min"),
         Rank_Weighted = rank(-Weighted_OutDegree, ties.method = "min"),
         Rank_Shift = Rank_Raw - Rank_Weighted) %>%
  arrange(Rank_Weighted)

cat("\nKnown FXS TFs — raw vs weighted (PANDA):\n")
print(as.data.frame(panda_strength %>% filter(TF %in% known_tf)))

write.csv(genie3_strength, "GENIE3_weighted_degree.csv", row.names = FALSE)
write.csv(panda_strength,  "PANDA_weighted_degree.csv",  row.names = FALSE)

known_rank_cmp <- bind_rows(
  genie3_strength %>% filter(TF %in% known_tf) %>%
    dplyr::select(TF, Rank_Raw, Rank_Weighted = Rank_MeanW) %>%
    mutate(Network = "GENIE3 (mean weight)"),
  panda_strength %>% filter(TF %in% known_tf) %>%
    dplyr::select(TF, Rank_Raw, Rank_Weighted) %>%
    mutate(Network = "PANDA (summed weight)"))

if (nrow(known_rank_cmp) > 0) {
  known_rank_cmp <- known_rank_cmp %>%
    pivot_longer(c(Rank_Raw, Rank_Weighted), names_to = "Metric",
                 values_to = "Rank") %>%
    mutate(Metric = recode(Metric, Rank_Raw = "Raw out-degree",
                           Rank_Weighted = "Weight-based"))

  ggplot(known_rank_cmp, aes(x = Metric, y = Rank, group = TF, color = TF)) +
    geom_line(linewidth = 1) + geom_point(size = 3) +
    geom_text_repel(aes(label = TF), size = 3.5, fontface = "bold",
                    max.overlaps = 20) +
    scale_y_reverse() + facet_wrap(~Network, scales = "free_y") +
    labs(title = "Known FXS TFs — Rank Under Raw vs Weight-Based Centrality",
         subtitle = "Higher on plot = more hub-like",
         x = "", y = "Rank among all TFs (reversed)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 9),
          legend.position = "none")
  ggsave("known_TF_raw_vs_weighted_rank.png", width = 10, height = 6, dpi = 300)
}

# ---- 18B: TARGET-LIST RESTRICTION BIAS ----

cat("\n--- 18B: Target-list restriction bias ---\n")

eligible_targets <- rownames(expr_subset)
all_measured     <- rownames(expr_normalized_sym)

restriction_bias <- data.frame(
  Metric = c("Genes in full expression matrix",
             "Genes eligible as GENIE3 targets",
             "Fraction of transcriptome eligible"),
  Value = c(length(all_measured), length(eligible_targets),
            round(length(eligible_targets) / length(all_measured), 4)))
print(restriction_bias)

motif_eligibility <- motif_edgelist %>%
  group_by(TF) %>%
  summarise(Motif_targets_total = n_distinct(Gene),
            Motif_targets_eligible = n_distinct(Gene[Gene %in% eligible_targets]),
            .groups = "drop") %>%
  mutate(Eligible_fraction = Motif_targets_eligible / Motif_targets_total)

cat("\nKnown FXS TFs — regulon eligibility in GENIE3 target space:\n")
print(as.data.frame(motif_eligibility %>% filter(TF %in% known_tf)))
cat("\nFOR DISCUSSION: a TF whose targets sit mostly outside the GENIE3\n")
cat("gene set cannot become a degree hub, regardless of importance.\n")

write.csv(motif_eligibility, "TF_target_eligibility_bias.csv", row.names = FALSE)
write.csv(restriction_bias, "GENIE3_restriction_bias_summary.csv", row.names = FALSE)

# ---- 18C: MODULE-LOCAL CENTRALITY ----
# PageRank is near-constant within these bipartite subgraphs (TFs have
# no in-edges), so its percentile is meaningless. Weighted out-strength
# is the correct local hub measure.

cat("\n--- 18C: Module-local centrality ---\n")

alpaca_modules <- alpaca_df %>% dplyr::select(Gene, Module) %>% distinct()
module_local_centrality <- list()

for (m in unique(alpaca_modules$Module)) {
  mod_genes <- alpaca_modules$Gene[alpaca_modules$Module == m]
  if (length(mod_genes) < 5) next
  sub_edges <- edge_list_panda %>%
    filter(TF %in% mod_genes, Target %in% mod_genes)
  if (nrow(sub_edges) < 3) next

  g_mod <- graph_from_data_frame(
    data.frame(from = sub_edges$TF, to = sub_edges$Target,
               weight = sub_edges$Score), directed = TRUE)

  cent <- data.frame(
    Gene = V(g_mod)$name, Module = m, Module_size = length(mod_genes),
    Local_Degree = degree(g_mod, mode = "all"),
    Local_OutDegree = degree(g_mod, mode = "out"),
    Local_Strength = strength(g_mod, mode = "out", weights = E(g_mod)$weight),
    Local_Between = betweenness(g_mod, directed = TRUE, normalized = TRUE),
    Local_PageRank = page_rank(g_mod, weights = E(g_mod)$weight)$vector,
    stringsAsFactors = FALSE)
  cent$Local_PageRank_pct <- rank(cent$Local_PageRank) / nrow(cent)
  cent$Local_Strength_pct <- rank(cent$Local_Strength) / nrow(cent)
  module_local_centrality[[as.character(m)]] <- cent
}

if (length(module_local_centrality) > 0) {
  module_cent_df <- do.call(rbind, module_local_centrality)
  cat("Modules with computable local centrality:",
      length(module_local_centrality), "\n")

  pr_var <- module_cent_df %>%
    group_by(Module) %>%
    summarise(n = n(),
              distinct_pagerank = n_distinct(round(Local_PageRank, 12)),
              distinct_strength = n_distinct(round(Local_Strength, 6)),
              .groups = "drop")
  cat("\nDistinct centrality values per module:\n"); print(as.data.frame(pr_var))
  cat("PageRank is near-constant within modules (bipartite, no TF\n")
  cat("in-edges). Weighted out-strength is used for ranking.\n")

  known_local <- module_cent_df %>%
    filter(Gene %in% known_tf) %>%
    arrange(desc(Local_Strength_pct))

  cat("\nKnown FXS TFs — module-local centrality (weighted strength):\n")
  if (nrow(known_local) > 0) {
    print(as.data.frame(known_local[, c("Gene","Module","Module_size",
                                        "Local_OutDegree","Local_Strength",
                                        "Local_Strength_pct")]))
    cat("\nIn top decile of their own module:",
        sum(known_local$Local_Strength_pct >= 0.9), "/", nrow(known_local), "\n")
    cat("Median percentile:",
        round(median(known_local$Local_Strength_pct), 3), "\n")
    cat("\nKEY RESULT: known FXS regulators may not be genome-wide degree\n")
    cat("hubs while still being locally central within ALPACA differential\n")
    cat("modules. Those are different claims; the second is the\n")
    cat("biologically meaningful one.\n")
  } else {
    cat("None of the known TFs appear in any module subnetwork.\n")
    cat("Report this directly — it is a real negative result.\n")
  }

  write.csv(module_cent_df, "ALPACA_module_local_centrality.csv", row.names = FALSE)

  if (nrow(known_local) > 0) {
    ggplot(known_local,
           aes(x = reorder(paste0(Gene, " (M", Module, ")"), Local_Strength_pct),
               y = Local_Strength_pct, fill = Local_Strength_pct >= 0.9)) +
      geom_bar(stat = "identity") + coord_flip() +
      geom_hline(yintercept = 0.9, linetype = "dashed") +
      scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#95A5A6"),
                        name = "Top decile") +
      labs(title = "Known FXS TFs — Centrality Within Their ALPACA Module",
           subtitle = "Percentile rank of weighted out-strength",
           x = "TF (module)", y = "Local weighted-strength percentile") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5))
    ggsave("known_TF_module_local_centrality.png", width = 8, height = 5, dpi = 300)
  }
}

cat("Step 18 done\n")

# ============================================================
# STEP 19 — IN SILICO VALIDATION
# ============================================================

cat("\n=== STEP 19: In silico validation ===\n")

top_novel_hubs <- head(setdiff(genie3_strength$TF, known_tf), 10)
validation_tfs <- unique(c(known_tf, top_novel_hubs,
                           head(panda_strength$TF, 10)))
validation_tfs <- intersect(validation_tfs,
                            unique(c(edge_list_genie3$TF, edge_list_panda$TF)))
cat("TFs to validate:", paste(validation_tfs, collapse = ", "), "\n")

chea3_query <- function(gene_set, query_name = "query") {
  gene_set <- unique(toupper(as.character(gene_set)))
  gene_set <- gene_set[!is.na(gene_set) & gene_set != ""]
  if (length(gene_set) < 5) return(NULL)
  resp <- tryCatch(
    POST("https://maayanlab.cloud/chea3/api/enrich/",
         body = toJSON(list(query_name = query_name, gene_set = gene_set),
                       auto_unbox = TRUE),
         encode = "raw", content_type_json(), timeout(120)),
    error = function(e) { cat("  request failed:", conditionMessage(e), "\n"); NULL })
  if (is.null(resp) || status_code(resp) != 200) return(NULL)
  tryCatch(fromJSON(content(resp, as = "text", encoding = "UTF-8")),
           error = function(e) NULL)
}

get_tf_targets <- function(tf, source = c("GENIE3", "PANDA"), top_n = 200) {
  source <- match.arg(source)
  if (source == "GENIE3") {
    edge_list_genie3 %>% filter(TF == tf) %>% arrange(desc(Weight)) %>%
      head(top_n) %>% pull(Target) %>% as.character()
  } else {
    edge_list_panda %>% filter(TF == tf) %>% arrange(desc(Score)) %>%
      head(top_n) %>% pull(Target) %>% as.character()
  }
}

chea3_results <- list()
for (tf in validation_tfs) {
  for (src in c("GENIE3", "PANDA")) {
    tgts <- get_tf_targets(tf, src)
    if (length(tgts) < 5) next
    cat("  ChEA3:", tf, "(", src, ",", length(tgts), "targets )\n")
    res_c <- chea3_query(tgts, paste0(tf, "_", src))
    if (is.null(res_c)) next

    lib <- if (!is.null(res_c[["Integrated--meanRank"]]))
      res_c[["Integrated--meanRank"]] else res_c[[1]]
    if (is.null(lib) || !"TF" %in% colnames(lib)) next
    lib$TF <- toupper(as.character(lib$TF))
    sr <- which(lib$TF == toupper(tf))
    self_rank <- if (length(sr) == 0) NA_integer_ else sr[1]

    enc <- res_c[["ENCODE--ChIP-seq"]]
    enc_rank <- NA_integer_
    if (!is.null(enc) && "TF" %in% colnames(enc)) {
      er <- which(toupper(as.character(enc$TF)) == toupper(tf))
      enc_rank <- if (length(er) == 0) NA_integer_ else er[1]
    }

    chea3_results[[paste0(tf, "_", src)]] <- data.frame(
      TF = tf, Source = src, N_targets = length(tgts),
      Self_rank_integrated = self_rank, N_TFs_returned = nrow(lib),
      Self_rank_ENCODE = enc_rank, Top_hit = lib$TF[1],
      stringsAsFactors = FALSE)
    Sys.sleep(1)
  }
}

if (length(chea3_results) > 0) {
  chea3_df <- do.call(rbind, chea3_results) %>% arrange(Self_rank_integrated)
  cat("\n--- ChEA3 integrated rank ---\n")
  print(as.data.frame(chea3_df))
  cat("\nTop 10 (integrated):",
      sum(chea3_df$Self_rank_integrated <= 10, na.rm = TRUE), "/",
      nrow(chea3_df), "\n")
  write.csv(chea3_df, "ChEA3_validation_results.csv", row.names = FALSE)

  # ENCODE is direct binding evidence — reported separately. The
  # integrated rank pools libraries where co-occurrence effects can
  # understate true recovery.
  cat("\n--- ENCODE ChIP-seq self-recovery (direct binding evidence) ---\n")
  encode_hits <- chea3_df %>% filter(!is.na(Self_rank_ENCODE)) %>%
    arrange(Self_rank_ENCODE) %>%
    dplyr::select(TF, Source, N_targets, Self_rank_ENCODE, Self_rank_integrated)

  if (nrow(encode_hits) > 0) {
    print(as.data.frame(encode_hits))
    cat("\nENCODE self-rank <= 10:",
        sum(encode_hits$Self_rank_ENCODE <= 10), "/", nrow(encode_hits), "\n")
    cat("ENCODE self-rank <= 25:",
        sum(encode_hits$Self_rank_ENCODE <= 25), "/", nrow(encode_hits), "\n")
    cat("TFs with ENCODE coverage:", nrow(encode_hits), "/", nrow(chea3_df), "\n")
    write.csv(encode_hits, "ChEA3_ENCODE_selfrecovery.csv", row.names = FALSE)

    ggplot(encode_hits,
           aes(x = reorder(paste0(TF, " (", Source, ")"), -Self_rank_ENCODE),
               y = Self_rank_ENCODE, fill = Self_rank_ENCODE <= 10)) +
      geom_bar(stat = "identity") + coord_flip() +
      geom_hline(yintercept = 10, linetype = "dashed") +
      scale_fill_manual(values = c("TRUE" = "#2ECC71", "FALSE" = "#95A5A6"),
                        name = "Rank <= 10") +
      labs(title = "ENCODE ChIP-seq Self-Recovery of Inferred Regulons",
           subtitle = "Rank of each TF in ENCODE enrichment for its own predicted targets",
           x = "TF (network)", y = "ENCODE library rank") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 8))
    ggsave("ChEA3_ENCODE_validation.png", width = 9, height = 6, dpi = 300)
    cat("\nReport ENCODE ranks as the validation claim.\n")
  }
} else {
  cat("No ChEA3 results — check network access to maayanlab.cloud\n")
}

# ---- 19B: TRRUST ----
# The universe includes TFs so that GENIE3 predictions (which are
# largely TF-TF here) are testable rather than silently dropped.

trrust_file <- "trrust_rawdata.human.tsv"
if (!file.exists(trrust_file)) {
  cat("\nDownloading TRRUST...\n")
  tryCatch(download.file(
    "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv",
    trrust_file, quiet = TRUE),
    error = function(e) cat("  download failed — fetch manually\n"))
}

if (file.exists(trrust_file)) {
  trrust <- read.table(trrust_file, sep = "\t", header = FALSE,
                       stringsAsFactors = FALSE, quote = "",
                       col.names = c("TF", "Target", "Direction", "PMID"))
  cat("TRRUST interactions:", nrow(trrust), "\n")

  universe_genes <- unique(as.character(
    c(edge_list_genie3$Target, edge_list_panda$Target,
      edge_list_genie3$TF, edge_list_panda$TF)))
  cat("Testing universe:", length(universe_genes), "genes\n")

  cat("\n--- TRRUST coverage diagnostic ---\n")
  for (tf in intersect(validation_tfs, unique(trrust$TF))) {
    g_pred <- unique(get_tf_targets(tf, "GENIE3", top_n = 500))
    p_pred <- unique(get_tf_targets(tf, "PANDA",  top_n = 500))
    curated <- unique(trrust$Target[trrust$TF == tf])
    cat(sprintf("%-8s curated=%3d | GENIE3 in-universe=%3d | PANDA in-universe=%3d\n",
                tf, length(curated),
                length(intersect(g_pred, universe_genes)),
                length(intersect(p_pred, universe_genes))))
  }

  trrust_validation <- lapply(validation_tfs, function(tf) {
    curated <- unique(trrust$Target[trrust$TF == tf])
    curated_in_universe <- intersect(curated, universe_genes)
    if (length(curated_in_universe) < 3) return(NULL)
    do.call(rbind, lapply(c("GENIE3", "PANDA"), function(src) {
      pred <- intersect(unique(get_tf_targets(tf, src, top_n = 500)),
                        universe_genes)
      if (length(pred) < 3) return(NULL)
      ov <- length(intersect(pred, curated_in_universe))
      expct <- length(pred) * length(curated_in_universe) / length(universe_genes)
      data.frame(TF = tf, Source = src, N_predicted = length(pred),
                 N_curated = length(curated_in_universe), N_overlap = ov,
                 Expected = expct, Fold_enrich = ov / expct,
                 pvalue = phyper(ov - 1, length(curated_in_universe),
                                 length(universe_genes) - length(curated_in_universe),
                                 length(pred), lower.tail = FALSE),
                 stringsAsFactors = FALSE)
    }))
  })

  trrust_df <- do.call(rbind, trrust_validation)
  if (!is.null(trrust_df) && nrow(trrust_df) > 0) {
    trrust_df$FDR <- p.adjust(trrust_df$pvalue, method = "BH")
    trrust_df <- trrust_df %>% arrange(pvalue)
    cat("\n--- TRRUST curated-target enrichment ---\n")
    print(as.data.frame(trrust_df))
    cat("\nSignificant (FDR<0.05):",
        sum(trrust_df$FDR < 0.05, na.rm = TRUE), "/", nrow(trrust_df), "\n")
    write.csv(trrust_df, "TRRUST_validation_results.csv", row.names = FALSE)

    ggplot(trrust_df, aes(x = reorder(TF, Fold_enrich), y = Fold_enrich,
                          fill = FDR < 0.05)) +
      geom_bar(stat = "identity", position = "dodge") + coord_flip() +
      facet_wrap(~Source) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      scale_fill_manual(values = c("TRUE" = "#2ECC71", "FALSE" = "#95A5A6"),
                        name = "FDR < 0.05") +
      labs(title = "Predicted Regulons vs TRRUST Curated Targets",
           subtitle = "Fold enrichment over chance",
           x = "Transcription Factor", y = "Fold enrichment") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5))
    ggsave("TRRUST_validation.png", width = 10, height = 6, dpi = 300)
  } else {
    cat("No TF had enough curated targets for testing.\n")
  }
} else {
  cat("TRRUST validation skipped.\n")
}

# ---- 19C: SUMMARY FOR THE MANUSCRIPT ----

cat("\n=== VALIDATION SUMMARY ===\n")
cat("In silico validation performed:\n")
cat("  - ChEA3 integrated TF-target library self-recovery\n")
cat("  - ENCODE ChIP-seq self-recovery (the primary claim)\n")
cat("  - TRRUST curated TF-target hypergeometric enrichment\n")

cat("\n--- METHODS WORDING: COHORTS ---\n")
cat("Primary: GSE146339, iPSC-derived cerebral organoids, 6 FXS vs 6\n")
cat("control, day 28 and day 105, timepoint as covariate.\n")
cat("Validation: GSE218676, patient-derived NPC / B-lymphocyte / iPSC,\n")
cat("5 FXS vs 5 control, cell type as covariate. Both cohorts quantified\n")
cat("with the same NCBI reprocessing pipeline (GRCh38.p13, Entrez), so\n")
cat("quantification differences cannot explain non-replication.\n")

cat("\n--- METHODS WORDING: METHYLATION LAYER ---\n")
cat("Methylation (GSE41273, peripheral blood) comes from an independent\n")
cat("cohort with no shared subjects. Promoter methylation is used as a\n")
cat("population-level prior on regulatory accessibility, not as matched\n")
cat("per-subject multi-omics. Describe the work as a multi-omic\n")
cat("integrative framework; do not claim subject-level integration.\n")

cat("\n--- LIMITATIONS TO STATE ---\n")
cat("1. Cohorts differ in tissue (cerebral organoid vs NPC/B-cell/iPSC).\n")
cat("   Genes replicating across both are robust to tissue context;\n")
cat("   genes failing to replicate may be tissue-specific rather than\n")
cat("   false positives.\n")
cat("2. Methylation is from a third tissue (blood) and a third cohort.\n")
cat("3. Validation cohort: NCBI reprocessed 20 of 26 submitted samples;\n")
cat("   FXS_426_NPC is absent. Two donors (FXS_421, NL_18) each\n")
cat("   contribute two cell types, so effective n is below the line\n")
cat("   count.\n")
cat("4. LIONESS has no covariate mechanism, so developmental timepoint\n")
cat("   is unadjusted in the sample-specific networks.\n")
cat("5. In silico validation shows inferred edges recover known binding\n")
cat("   above chance. It does not establish causality in FXS and is not\n")
cat("   equivalent to experimental confirmation.\n")

save.image("checkpoint_step19.RData")
cat("\n=== PIPELINE COMPLETE ===\n")
