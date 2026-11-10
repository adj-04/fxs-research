# ============================================================
# FXS PIPELINE v5 — PART 2 of 3 : STEPS 9–16
# Networks, hubs, LIONESS, ALPACA, enrichment, Cytoscape,
# validation cohort, cross-cohort concordance
#
# Requires: checkpoint_step8.RData (from part 1)
# Ends by saving checkpoint_step16.RData
# Then: restart R and run part3_steps17-19.R
#
# v5 changes: memory-safe LIONESS matrix; per-edge limma test on LIONESS
# networks (~ timepoint + condition); GO/KEGG universe fixed for relaxed
# runs; Cytoscape-free XGMML export; donor-aware limma-voom sensitivity
# for the validation cohort; shrunken-LFC cross-cohort correlation.
# ============================================================

library(DESeq2)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(scales)
library(Biobase)
library(org.Hs.eg.db)
library(clusterProfiler)
library(pandaR)
library(netZooR)
library(igraph)
if (toupper(Sys.getenv("FXS_RUN_CYTOSCAPE", "FALSE")) == "TRUE") library(RCy3)
library(limma)
library(edgeR)
library(ggvenn)

# Paths/config, then the checkpoint, then config again so that settings
# from the environment (not the values saved inside the checkpoint) win.
source(file.path(Sys.getenv("FXS_REPO", "."), "scripts", "00_config.R"))
load(ckpt("checkpoint_step8.RData"))
source(file.path(Sys.getenv("FXS_REPO", "."), "scripts", "00_config.R"))
cat("Checkpoint loaded. Objects:", length(ls()), "\n")
log_session("part2")

# Cytoscape cannot run on Colab. v5 writes XGMML/GraphML directly from R;
# set FXS_RUN_CYTOSCAPE=TRUE locally to also drive a running Cytoscape.
RUN_CYTOSCAPE <- env_bool("FXS_RUN_CYTOSCAPE", "FALSE")
VALID_EXCLUDE_PREMUTATION <- TRUE
VALID_MODE                <- "pooled"
set.seed(42)

# Steps 9-12 (networks, LIONESS, ALPACA) are the slow part. They are
# checkpointed, so a failure later in Part 2 does not repeat them.
# Set FXS_REDO_STEP12=TRUE to force a rerun.
resume12 <- file.exists(ckpt("checkpoint_step12.RData")) &&
  !env_bool("FXS_REDO_STEP12", "FALSE") &&
  file.mtime(ckpt("checkpoint_step12.RData")) > file.mtime(ckpt("checkpoint_step8.RData"))
if (resume12) {
  cat("Resuming from checkpoint_step12.RData (Steps 9-12 already done)\n")
  load(ckpt("checkpoint_step12.RData"))
  source(file.path(Sys.getenv("FXS_REPO", "."), "scripts", "00_config.R"))
} else {

# ============================================================
# STEP 9 — NETWORK CONSTRUCTION WITH IGRAPH
# ============================================================

cat("\nBuilding GENIE3 network...\n")
genie3_threshold <- quantile(edge_list_genie3$Weight, GENIE3_EDGE_QUANTILE)
genie3_filtered  <- edge_list_genie3[edge_list_genie3$Weight >= genie3_threshold, ]
cat("GENIE3 threshold (", 100 * GENIE3_EDGE_QUANTILE, "th pct):", signif(genie3_threshold, 4),
    "| edges:", nrow(genie3_filtered), "\n")

g_genie3 <- graph_from_data_frame(
  data.frame(from = genie3_filtered$TF, to = genie3_filtered$Target,
             weight = genie3_filtered$Weight), directed = TRUE)
cat("GENIE3 network — Nodes:", vcount(g_genie3),
    "Edges:", ecount(g_genie3), "\n")

cat("\nBuilding PANDA network...\n")
panda_threshold <- quantile(edge_list_panda$Score, 0.90)
panda_filtered  <- edge_list_panda[edge_list_panda$Score >= panda_threshold, ]
cat("PANDA threshold (90th pct):", signif(panda_threshold, 4),
    "| edges:", nrow(panda_filtered), "\n")

g_panda <- graph_from_data_frame(
  data.frame(from = panda_filtered$TF, to = panda_filtered$Target,
             score = panda_filtered$Score), directed = TRUE)
cat("PANDA network — Nodes:", vcount(g_panda),
    "Edges:", ecount(g_panda), "\n")

cat("\n--- Network statistics ---\n")
cat("GENIE3 density:", signif(edge_density(g_genie3), 4),
    "| components:", igraph::components(g_genie3)$no, "\n")
cat("PANDA density:", signif(edge_density(g_panda), 4),
    "| components:", igraph::components(g_panda)$no, "\n")

V(g_genie3)$type <- ifelse(V(g_genie3)$name %in% tf_in_data, "TF", "Target")
V(g_panda)$type  <- ifelse(V(g_panda)$name  %in% tf_in_data, "TF", "Target")
cat("GENIE3 TF/Target nodes:", sum(V(g_genie3)$type == "TF"), "/",
    sum(V(g_genie3)$type == "Target"), "\n")
cat("PANDA TF/Target nodes:", sum(V(g_panda)$type == "TF"), "/",
    sum(V(g_panda)$type == "Target"), "\n")

write_graph(g_genie3, file.path(RESULTS_GENIE3, "GENIE3_network.graphml"), format = "graphml")
write_graph(g_panda,  file.path(RESULTS_PANDA, "PANDA_network.graphml"),  format = "graphml")

for (nm in c("genie3", "panda")) {
  g   <- if (nm == "genie3") g_genie3 else g_panda
  col <- if (nm == "genie3") "#3498DB" else "#2ECC71"
  fig_dir <- if (nm == "genie3") FIG_GENIE3 else FIG_PANDA
  V(g)$degree <- degree(g, mode = "all")
  set.seed(42); lay <- layout_with_fr(g)
  png(file.path(fig_dir, paste0(toupper(nm), "_network_plot.png")),
      width = 2000, height = 2000, res = 200)
  plot(g, layout = lay,
       vertex.size = scales::rescale(V(g)$degree, to = c(1, 8)),
       vertex.color = ifelse(V(g)$type == "TF", "#E74C3C", col),
       vertex.frame.color = NA,
       vertex.label = ifelse(V(g)$degree > quantile(V(g)$degree, 0.95),
                             V(g)$name, NA),
       vertex.label.cex = 0.5, vertex.label.color = "black",
       edge.arrow.size = 0.1,
       edge.color = adjustcolor("grey50", alpha.f = 0.3), edge.width = 0.3,
       main = paste(toupper(nm), "Gene Regulatory Network — FXS"))
  legend("bottomleft", legend = c("TF", "Target Gene"),
         col = c("#E74C3C", col), pch = 19, pt.cex = 1.5, bty = "n")
  dev.off()
}
cat("Network plots saved\n")

# ============================================================
# STEP 10 — HUB GENE IDENTIFICATION
# ============================================================

cat("\n--- Step 10: Hub gene identification ---\n")

genie3_centrality <- data.frame(
  Gene = V(g_genie3)$name, Type = V(g_genie3)$type,
  Degree      = degree(g_genie3, mode = "all"),
  InDegree    = degree(g_genie3, mode = "in"),
  OutDegree   = degree(g_genie3, mode = "out"),
  Betweenness = betweenness(g_genie3, directed = TRUE, normalized = TRUE),
  PageRank    = page_rank(g_genie3)$vector,
  OutStrength = strength(g_genie3, mode = "out", weights = E(g_genie3)$weight))
# v5: hubs = TFs ranked by weighted out-degree (as in the paper, Fig. 7).
# PageRank on a directed TF->target graph ranks receivers, not regulators.
genie3_centrality <- genie3_centrality[order(genie3_centrality$OutStrength,
                                             decreasing = TRUE), ]
cat("\nTop 20 hub TFs — GENIE3 (by weighted out-degree):\n")
print(head(genie3_centrality, 20))
write.csv(genie3_centrality, file.path(RESULTS_GENIE3, "GENIE3_centrality.csv"), row.names = FALSE)

panda_centrality <- data.frame(
  Gene = V(g_panda)$name, Type = V(g_panda)$type,
  Degree      = degree(g_panda, mode = "all"),
  InDegree    = degree(g_panda, mode = "in"),
  OutDegree   = degree(g_panda, mode = "out"),
  Betweenness = betweenness(g_panda, directed = TRUE, normalized = TRUE),
  PageRank    = page_rank(g_panda)$vector)
panda_centrality <- panda_centrality[order(panda_centrality$OutDegree,
                                           decreasing = TRUE), ]
cat("\nTop 20 hub genes — PANDA (by OutDegree):\n")
print(head(panda_centrality, 20))
write.csv(panda_centrality, file.path(RESULTS_PANDA, "PANDA_centrality.csv"), row.names = FALSE)

cat("\nKnown FXS TFs in GENIE3 top 50:\n")
print(intersect(known_tf, head(genie3_centrality$Gene, 50)))
cat("\nKnown FXS TFs in PANDA top 50:\n")
print(intersect(known_tf, head(panda_centrality$Gene, 50)))

top20_genie3 <- head(genie3_centrality, 20)
top20_genie3$Gene <- factor(top20_genie3$Gene,
                            levels = top20_genie3$Gene[order(top20_genie3$OutStrength)])
ggplot(top20_genie3, aes(x = Gene, y = OutStrength, fill = Type)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = c("TF" = "#E74C3C", "Target" = "#3498DB")) +
  labs(title = "Top 20 Hub TFs - GENIE3 Network (weighted out-degree)",
       x = "TF", y = "Weighted out-degree") + theme_bw()
ggsave(file.path(FIG_GENIE3, "GENIE3_hub_genes.png"), width = 8, height = 6, dpi = 300)

top20_panda <- head(panda_centrality, 20)
top20_panda$Gene <- factor(top20_panda$Gene,
                           levels = top20_panda$Gene[order(top20_panda$OutDegree)])
ggplot(top20_panda, aes(x = Gene, y = OutDegree, fill = Type)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = c("TF" = "#E74C3C", "Target" = "#2ECC71")) +
  labs(title = "Top 20 Hub TFs — PANDA Network (by OutDegree)",
       x = "Gene", y = "OutDegree (targets regulated)") + theme_bw()
ggsave(file.path(FIG_PANDA, "PANDA_hub_genes.png"), width = 8, height = 6, dpi = 300)

genie3_known <- genie3_centrality[genie3_centrality$Gene %in% known_tf,
                                  c("Gene","OutDegree","Betweenness","PageRank")]
panda_known  <- panda_centrality[panda_centrality$Gene %in% known_tf,
                                 c("Gene","OutDegree","Betweenness","PageRank")]
genie3_known$Network <- "GENIE3"    # safe on zero rows; keeps columns aligned
panda_known$Network  <- "PANDA"
combined_known <- rbind(genie3_known, panda_known)

if (nrow(combined_known) > 0) {
  ggplot(combined_known, aes(x = Gene, y = OutDegree, fill = Network)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    scale_fill_manual(values = c("GENIE3" = "#E74C3C", "PANDA" = "#3498DB")) +
    labs(title = "Known FXS Regulators — OutDegree in GENIE3 vs PANDA",
         subtitle = "Higher OutDegree = more target genes regulated",
         x = "Transcription Factor", y = "OutDegree") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  ggsave(file.path(FIG_OTHER, "FXS_TF_OutDegree_comparison.png"), width = 9, height = 6, dpi = 300)
}

cat("Step 10 done\n")

# ============================================================
# STEP 11 — LIONESS SAMPLE-SPECIFIC NETWORKS
# ============================================================

cat("\n=== Applying netZooR/pandaR class() patches ===\n")

fix_class_check <- function(fn) {
  fix_node <- function(e) {
    if (is.call(e)) {
      if (identical(e[[1]], as.name("if"))) {
        cond <- e[[2]]
        if (is.call(cond) && identical(cond[[1]], as.name("==")) &&
            is.call(cond[[2]]) && identical(cond[[2]][[1]], as.name("class")) &&
            identical(cond[[3]], "ExpressionSet")) {
          e[[2]] <- call("inherits", cond[[2]][[2]], "ExpressionSet")
        }
      }
      for (i in seq_along(e)) e[[i]] <- fix_node(e[[i]])
    }
    e
  }
  body(fn) <- fix_node(body(fn)); fn
}

patch_namespace <- function(ns, ns_name) {
  for (fn_name in ls(ns)) {
    obj <- get(fn_name, envir = ns)
    if (is.function(obj)) {
      b <- tryCatch(deparse(body(obj)), error = function(e) "")
      if (any(grepl('class\\(.*\\) == "ExpressionSet"', b))) {
        fixed <- fix_class_check(obj)
        unlockBinding(fn_name, ns); assign(fn_name, fixed, envir = ns)
        lockBinding(fn_name, ns)
        cat("  Patched:", ns_name, "::", fn_name, "\n")
      }
    }
  }
}

patch_namespace(asNamespace("pandaR"),  "pandaR")
patch_namespace(asNamespace("netZooR"), "netZooR")

expr_lioness <- as.matrix(expr_final)
storage.mode(expr_lioness) <- "double"
cat("\nLIONESS input:", dim(expr_lioness), "\n")
if (ncol(expr_lioness) < 2) stop("Need at least 2 samples for LIONESS")

expr_es <- new("ExpressionSet", exprs = expr_lioness)
motif_lioness <- as.data.frame(motif_edgelist)
colnames(motif_lioness) <- c("TF", "Gene", "Score")
ppi_lioness <- as.data.frame(ppi_edgelist)
colnames(ppi_lioness) <- c("TF1", "TF2", "Score")

cat("Running LIONESS —", ncol(expr_lioness), "networks...\n")
lioness_result <- netZooR::lioness(
  expr = expr_es, motif = motif_lioness, ppi = ppi_lioness,
  network.inference.method = "panda", ncores = 1, mode = "union")
cat("LIONESS complete:", length(lioness_result), "networks\n")

sample_names_lioness <- colnames(expr_final)

# v5: build the edge x sample matrix column by column. v4 built a long
# data frame (edges x samples rows of pasted strings) and dcast it,
# which needs several GB once the target universe grows.
extract_lioness <- function(x) {
  if (is.data.frame(x)) {
    list(id = paste(as.character(x$TF), as.character(x$Gene), sep = "__"),
         score = as.numeric(as.character(x$Score)))
  } else {
    x <- as.matrix(x)
    list(id = paste(rep(rownames(x), times = ncol(x)),
                    rep(colnames(x), each = nrow(x)), sep = "__"),
         score = as.numeric(x))
  }
}
first <- extract_lioness(lioness_result[[1]])
lioness_mat <- matrix(NA_real_, nrow = length(first$id),
                      ncol = length(lioness_result),
                      dimnames = list(first$id, sample_names_lioness))
lioness_mat[, 1] <- first$score
for (i in seq_along(lioness_result)[-1]) {
  e <- extract_lioness(lioness_result[[i]])
  lioness_mat[, i] <- if (identical(e$id, first$id)) e$score else
    e$score[match(first$id, e$id)]
  lioness_result[i] <- list(NULL)   # free each network once copied
  rm(e)
}
rm(lioness_result, first); gc(verbose = FALSE)
lioness_mat <- lioness_mat[stats::complete.cases(lioness_mat), , drop = FALSE]
lioness_df  <- lioness_mat          # name kept for compatibility with part 3
cat("LIONESS edge matrix:", dim(lioness_mat),
    sprintf("(%.0f MB)\n", as.numeric(object.size(lioness_mat)) / 1e6))

fxs_cols     <- intersect(names(condition)[condition == "FXS"],
                          colnames(lioness_mat))
control_cols <- intersect(names(condition)[condition == "Control"],
                          colnames(lioness_mat))
cat("FXS:", length(fxs_cols), "| Control:", length(control_cols), "\n")
cat("FXS:", paste(fxs_cols, collapse = ", "), "\n")
cat("Control:", paste(control_cols, collapse = ", "), "\n")
if (length(fxs_cols) == 0 || length(control_cols) == 0)
  stop("Sample name mismatch — check names(condition) vs colnames(lioness_mat)")

fxs_mean     <- rowMeans(lioness_mat[, fxs_cols,     drop = FALSE])
control_mean <- rowMeans(lioness_mat[, control_cols, drop = FALSE])

# ---- v5: per-edge differential test with the timepoint covariate ----
# LIONESS itself cannot take covariates, but its output can: each edge's
# 12 sample-specific weights are modelled as ~ timepoint + condition
# with limma (the approach used in the LIONESS/netZoo papers). The input
# expression was already timepoint-adjusted (Step 6); keeping timepoint
# in the edge model removes any residual non-linear stage effect.
edge_design <- model.matrix(~ timepoint[colnames(lioness_mat)] +
                              condition[colnames(lioness_mat)])
colnames(edge_design) <- c("Intercept", "timepointD105", "conditionFXS")
edge_fit <- eBayes(lmFit(lioness_mat, edge_design))
edge_tt  <- topTable(edge_fit, coef = "conditionFXS", number = Inf,
                     sort.by = "none")
rm(edge_fit); gc(verbose = FALSE)

diff_df <- data.frame(Edge = rownames(lioness_mat), FXS = fxs_mean,
                      Control = control_mean, Diff = fxs_mean - control_mean,
                      limma_logFC = edge_tt$logFC, limma_t = edge_tt$t,
                      limma_p = edge_tt$P.Value, limma_FDR = edge_tt$adj.P.Val,
                      stringsAsFactors = FALSE)
rm(edge_tt)
diff_df$TF     <- sub("__.*", "", diff_df$Edge)
diff_df$Target <- sub(".*__", "", diff_df$Edge)
diff_df <- diff_df[order(diff_df$limma_p), ]

cat("Differential edges:", nrow(diff_df),
    "| FXS >:", sum(diff_df$Diff > 0), "| FXS <:", sum(diff_df$Diff < 0), "\n")
cat("Edges, condition effect adjusted for timepoint — nominal p<0.05:",
    sum(diff_df$limma_p < 0.05), "| FDR<0.05:", sum(diff_df$limma_FDR < 0.05),
    "| FDR<0.10:", sum(diff_df$limma_FDR < 0.10), "\n")
write.csv(diff_df, file.path(RESULTS_LIONESS, "LIONESS_diff_edges.csv"), row.names = FALSE)
write.csv(head(diff_df, 5000),
          file.path(RESULTS_LIONESS, "LIONESS_diff_edges_top5000_limma.csv"), row.names = FALSE)

# p-value histogram: uniform = no signal; spike near 0 = signal
png(file.path(FIG_LIONESS, "LIONESS_limma_pvalue_hist.png"), width = 1600, height = 1100, res = 200)
hist(diff_df$limma_p, breaks = 50, col = "#3498DB", border = "white",
     main = "LIONESS per-edge test (~ timepoint + condition)",
     xlab = "Nominal p-value, condition effect")
dev.off()

plot_idx <- if (nrow(diff_df) > 200000) sample(nrow(diff_df), 200000) else seq_len(nrow(diff_df))
ggplot(diff_df[plot_idx, ], aes(x = Control, y = FXS)) +
  geom_point(alpha = 0.3, size = 0.5, color = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "LIONESS Edge Score Comparison — FXS vs Control",
       subtitle = "Points above the dashed line are stronger in FXS",
       x = "Control network score", y = "FXS network score") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
ggsave(file.path(FIG_LIONESS, "LIONESS_FXS_vs_Control_scatter.png"), width = 8, height = 6, dpi = 300)

top_edges <- head(diff_df[order(abs(diff_df$Diff), decreasing = TRUE), ], 20)
top_edges$Label <- paste0(top_edges$TF, " -> ", top_edges$Target)
ggplot(top_edges, aes(x = reorder(Label, Diff), y = Diff, fill = Diff > 0)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#3498DB"),
                    labels = c("TRUE" = "Increased in FXS",
                               "FALSE" = "Decreased in FXS"),
                    name = "Direction") +
  labs(title = "Top 20 Differential TF-Target Edges (FXS vs Control)",
       subtitle = "Ranked by absolute differential score",
       x = "Edge", y = "Differential score (FXS - Control)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
ggsave(file.path(FIG_LIONESS, "LIONESS_top_edges_barplot.png"), width = 9, height = 6, dpi = 300)

known_diff <- diff_df[diff_df$TF %in% known_tf, ]
if (nrow(known_diff) > 0) {
  known_summary <- known_diff %>%
    group_by(TF) %>%
    summarise(Mean_Diff = mean(Diff), N_edges = n(), .groups = "drop")
  print(known_summary)
  ggplot(known_summary, aes(x = reorder(TF, Mean_Diff), y = Mean_Diff,
                            fill = Mean_Diff > 0)) +
    geom_bar(stat = "identity") + coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#3498DB"),
                      labels = c("TRUE" = "Hyperactive in FXS",
                                 "FALSE" = "Hypoactive in FXS"),
                      name = "Activity shift") +
    labs(title = "Known FXS TF Activity Shift (LIONESS)",
         subtitle = "Mean differential edge score across all targets",
         x = "Transcription Factor", y = "Mean differential score") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.y = element_text(size = 12, face = "bold"))
  ggsave(file.path(FIG_LIONESS, "LIONESS_known_TF_activity.png"), width = 7, height = 5, dpi = 300)
}

tfs   <- sub("__.*", "", rownames(lioness_mat))
genes <- sub(".*__", "", rownames(lioness_mat))

fxs_net     <- data.frame(TF = tfs, Gene = genes, Score = fxs_mean,
                          stringsAsFactors = FALSE)
control_net <- data.frame(TF = tfs, Gene = genes, Score = control_mean,
                          stringsAsFactors = FALSE)
cat("Group networks:", nrow(fxs_net), "edges each\n")
cat("TFs:", length(unique(fxs_net$TF)),
    "| Targets:", length(unique(fxs_net$Gene)), "\n")
if (length(unique(fxs_net$TF)) >= length(unique(fxs_net$Gene)))
  cat("WARNING: TFs >= targets — condor/ALPACA will fail. Raise N_TARGET_POOL.\n")

write.csv(fxs_net,     file.path(RESULTS_ALPACA, "FXS_group_network.csv"),     row.names = FALSE)
write.csv(control_net, file.path(RESULTS_OTHER, "Control_group_network.csv"), row.names = FALSE)
cat("Step 11 done\n")

# ============================================================
# STEP 12 — ALPACA DIFFERENTIAL MODULE DETECTION
# ============================================================

cat("\nRunning ALPACA...\n")
# v5.1 FIX: alpaca() treats the 3rd column of the merged table (panda.net1)
# as the CONTROL network (communities are detected in it), so control goes
# first. v5 passed FXS first. Note the module scores written by alpaca()
# sum to 2 within every module; module significance is tested in part 5.
alpaca_result <- netZooR::pandaToAlpaca(
  panda.net1 = control_net, panda.net2 = fxs_net,
  file.stem = file.path(RESULTS_ALPACA, "FXS_ALPACA"), verbose = TRUE)

memb_file <- read.table(file.path(RESULTS_ALPACA, "FXS_ALPACA_ALPACA_final_memb.txt"), header = FALSE,
                        sep = "\t", stringsAsFactors = FALSE,
                        col.names = c("Node", "Module"))
scores_file <- read.table(file.path(RESULTS_ALPACA, "FXS_ALPACA_ALPACA_scores.txt"), header = FALSE,
                          sep = "\t", stringsAsFactors = FALSE,
                          col.names = c("Node", "Score"))

alpaca_df <- merge(memb_file, scores_file, by = "Node", all.x = TRUE)
alpaca_df$Score[is.na(alpaca_df$Score)] <- 0
alpaca_df$Gene <- gsub("_A$|_B$", "", alpaca_df$Node)
alpaca_df$Type <- ifelse(grepl("_A$", alpaca_df$Node), "TF", "Target")
alpaca_df <- alpaca_df[order(alpaca_df$Score, decreasing = TRUE), ]

cat("ALPACA nodes:", nrow(alpaca_df),
    "| modules:", length(unique(alpaca_df$Module)), "\n")
cat("\nTop 20 ALPACA nodes:\n")
print(head(alpaca_df[, c("Gene", "Module", "Score", "Type")], 20))
cat("\nKnown FXS TFs in ALPACA:\n")
print(alpaca_df[alpaca_df$Gene %in% known_tf,
                c("Gene", "Module", "Score", "Type")])
write.csv(alpaca_df, file.path(RESULTS_ALPACA, "ALPACA_results.csv"), row.names = FALSE)

module_counts <- alpaca_df %>% count(Module, name = "Count") %>%
  arrange(desc(Count))
ggplot(head(module_counts, 15),
       aes(x = reorder(Module, Count), y = Count, fill = Count)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_gradient(low = "#3498DB", high = "#E74C3C") +
  labs(title = "ALPACA Module Sizes — FXS vs Control",
       x = "Module", y = "Nodes") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(FIG_ALPACA, "ALPACA_module_sizes.png"), width = 8, height = 6, dpi = 300)

known_alpaca <- alpaca_df[alpaca_df$Gene %in% known_tf, ]
if (nrow(known_alpaca) > 0) {
  ggplot(known_alpaca, aes(x = reorder(Gene, Score), y = Score,
                           fill = factor(Module))) +
    geom_bar(stat = "identity") + coord_flip() +
    scale_fill_brewer(palette = "Set1", name = "Module") +
    labs(title = "Known FXS TFs — ALPACA Module Assignment",
         x = "Transcription Factor", y = "ALPACA Score") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(FIG_ALPACA, "ALPACA_FXS_TF_modules.png"), width = 8, height = 5, dpi = 300)
}

ggplot(head(alpaca_df, 20), aes(x = reorder(Gene, Score), y = Score, fill = Type)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = c("TF" = "#E74C3C", "Target" = "#3498DB")) +
  labs(title = "Top 20 ALPACA Nodes — FXS Differential Modules",
       x = "Gene", y = "ALPACA Score") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(FIG_ALPACA, "ALPACA_top_nodes.png"), width = 8, height = 6, dpi = 300)

cat("Step 12 done\n")
save.image(ckpt("checkpoint_step12.RData"))
cat("checkpoint_step12.RData saved\n")
}   # end of Steps 9-12 block

# ============================================================
# STEP 13 — FUNCTIONAL ENRICHMENT
# ============================================================

alpaca_df <- alpaca_df[, !duplicated(colnames(alpaca_df))]

module_sizes <- alpaca_df %>% group_by(Module) %>%
  summarise(n = n(), .groups = "drop") %>% filter(n >= 5)
cat("Modules with >= 5 genes:", nrow(module_sizes), "\n")

gene_map <- bitr(unique(alpaca_df$Gene), fromType = "SYMBOL",
                 toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat("Genes mapped to Entrez:", nrow(gene_map), "/",
    length(unique(alpaca_df$Gene)), "\n")

alpaca_df <- alpaca_df[, !colnames(alpaca_df) %in% "ENTREZID"]
alpaca_df <- merge(alpaca_df, gene_map, by.x = "Gene", by.y = "SYMBOL",
                   all.x = TRUE)

cat("\nRunning per-module GO enrichment...\n")
module_go_results <- compareCluster(
  ENTREZID ~ Module,
  data = alpaca_df %>%
    filter(Module %in% module_sizes$Module, !is.na(ENTREZID)) %>%
    dplyr::select(Module, ENTREZID),
  fun = "enrichGO", OrgDb = org.Hs.eg.db, ont = "BP",
  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05,
  readable = TRUE)

go_df <- as.data.frame(module_go_results)
write.csv(go_df, file.path(RESULTS_ALPACA, "ALPACA_module_GO_enrichment.csv"), row.names = FALSE)

if (nrow(go_df) > 0) {
  dotplot(module_go_results, showCategory = 5, includeAll = FALSE) +
    ggtitle("GO BP Enrichment — ALPACA Modules") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(FIG_ALPACA, "ALPACA_module_GO_dotplot.png"), width = 12, height = 8, dpi = 300)
}

top_hits <- alpaca_df %>% filter(Score > 0.05, !is.na(ENTREZID)) %>%
  pull(ENTREZID) %>% unique()
background <- alpaca_df %>% filter(!is.na(ENTREZID)) %>%
  pull(ENTREZID) %>% unique()
cat("Top hit genes (Score > 0.05):", length(top_hits), "\n")

top_go_strict <- enrichGO(gene = top_hits, universe = background,
                          OrgDb = org.Hs.eg.db, ont = "BP",
                          pAdjustMethod = "BH", pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05, readable = TRUE)
top_kegg_strict <- tryCatch(
  enrichKEGG(gene = top_hits, universe = background, organism = "hsa",
             pAdjustMethod = "BH", pvalueCutoff = 0.05),
  error = function(e) { cat("KEGG strict failed:", conditionMessage(e), "\n"); NULL })

if (!is.null(top_go_strict)) {
  cat("GO terms (strict):", nrow(as.data.frame(top_go_strict)), "\n")
  if (nrow(as.data.frame(top_go_strict)) > 0)
    print(head(as.data.frame(top_go_strict)[, c("Description","GeneRatio",
                                                "pvalue","p.adjust")], 20))
}

# v5: the relaxed runs now use the same universe as the strict runs. v4
# omitted it, so relaxed terms were tested against the whole genome and
# partly reflected which genes could enter the network at all.
top_go_relaxed <- enrichGO(gene = top_hits, universe = background,
                           OrgDb = org.Hs.eg.db, ont = "BP",
                           pAdjustMethod = "BH", pvalueCutoff = 0.1,
                           qvalueCutoff = 0.2, readable = TRUE)
top_kegg_relaxed <- tryCatch(
  enrichKEGG(gene = top_hits, universe = background, organism = "hsa",
             pAdjustMethod = "BH", pvalueCutoff = 0.1),
  error = function(e) NULL)

if (!is.null(top_go_relaxed))
  cat("GO terms (relaxed):", nrow(as.data.frame(top_go_relaxed)), "\n")

if (!is.null(top_go_strict))
  write.csv(as.data.frame(top_go_strict),  file.path(RESULTS_ALPACA, "ALPACA_topHits_GO_strict.csv"),  row.names = FALSE)
if (!is.null(top_go_relaxed))
  write.csv(as.data.frame(top_go_relaxed), file.path(RESULTS_ALPACA, "ALPACA_topHits_GO_relaxed.csv"), row.names = FALSE)
if (!is.null(top_kegg_strict))
  write.csv(as.data.frame(top_kegg_strict),  file.path(RESULTS_ALPACA, "ALPACA_topHits_KEGG_strict.csv"),  row.names = FALSE)
if (!is.null(top_kegg_relaxed))
  write.csv(as.data.frame(top_kegg_relaxed), file.path(RESULTS_ALPACA, "ALPACA_topHits_KEGG_relaxed.csv"), row.names = FALSE)

best_go <- top_go_relaxed
if (!is.null(top_go_strict) && nrow(as.data.frame(top_go_strict)) >= 3) {
  best_go <- top_go_strict
}

if (!is.null(best_go) && nrow(as.data.frame(best_go)) > 0) {
  dotplot(best_go, showCategory = 20) +
    ggtitle("GO BP — Top ALPACA Scoring Genes") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(FIG_ALPACA, "ALPACA_topHits_GO_dotplot.png"), width = 10, height = 8, dpi = 300)
}

best_kegg <- NULL
if (!is.null(top_kegg_strict) &&
    nrow(as.data.frame(top_kegg_strict)) >= 3) {
  best_kegg <- top_kegg_strict
} else if (!is.null(top_kegg_relaxed) &&
           nrow(as.data.frame(top_kegg_relaxed)) > 0) {
  best_kegg <- top_kegg_relaxed
}

if (!is.null(best_kegg)) {
  dotplot(best_kegg, showCategory = 20) +
    ggtitle("KEGG Pathways — Top ALPACA Scoring Genes") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(FIG_ALPACA, "ALPACA_topHits_KEGG_dotplot.png"), width = 10, height = 8, dpi = 300)
}
cat("Step 13 done\n")

# ============================================================
# STEP 14 — CYTOSCAPE EXPORT
# ============================================================

if (RUN_CYTOSCAPE) {
  
  top_module_ids <- alpaca_df %>%
    group_by(Module) %>%
    summarise(mean_score = mean(Score), n = n(), .groups = "drop") %>%
    filter(n >= 5) %>% arrange(desc(mean_score)) %>% head(3) %>% pull(Module)
  
  top_nodes <- alpaca_df %>%
    filter(Module %in% top_module_ids) %>%
    dplyr::select(Gene, Module, Score, Type) %>%
    distinct()
  cat("Top modules:", top_module_ids, "| nodes:", nrow(top_nodes), "\n")
  
  edge_table <- edge_list_panda %>%
    filter(TF %in% top_nodes$Gene, Target %in% top_nodes$Gene) %>%
    rename(Gene = Target)
  cat("Edges in subnetwork:", nrow(edge_table), "\n")
  
  cytoscapePing()
  
  node_table <- top_nodes %>%
    rename(id = Gene) %>%
    mutate(label     = id,
           Module    = as.character(Module),
           Score     = round(Score, 4),
           NodeColor = ifelse(Type == "TF", "#E74C3C", "#3498DB"))
  
  createNetworkFromDataFrames(
    nodes = node_table,
    edges = edge_table %>% rename(source = TF, target = Gene),
    title = "FXS_ALPACA_TopModules", collection = "FXS_PANDA_ALPACA")
  cat("Network created\n")
  
  cat("Node table columns in Cytoscape:\n")
  print(head(getTableColumns("node"), 3))
  
  style_name <- "FXS_ALPACA_Style"
  if (style_name %in% getVisualStyleNames()) deleteVisualStyle(style_name)
  createVisualStyle(style_name)
  
  setNodeColorMapping(table.column = "Type",
                      table.column.values = c("TF", "Target"),
                      colors = c("#E74C3C", "#3498DB"),
                      mapping.type = "d", style.name = style_name)
  setNodeSizeMapping(table.column = "Score",
                     table.column.values = c(0, 0.1, 0.3, 1.0),
                     sizes = c(35, 50, 70, 100),
                     mapping.type = "c", style.name = style_name)
  setNodeLabelMapping(table.column = "label", style.name = style_name)
  setNodeShapeMapping(table.column = "Type",
                      table.column.values = c("TF", "Target"),
                      shapes = c("ELLIPSE", "ROUND_RECTANGLE"),
                      style.name = style_name)
  
  setNodeBorderColorDefault("#FFFFFF", style.name = style_name)
  setNodeBorderWidthDefault(2,  style.name = style_name)
  setNodeFontSizeDefault(20, style.name = style_name)
  setNodeLabelColorDefault("#000000", style.name = style_name)
  setNodeLabelOpacityDefault(255, style.name = style_name)
  setEdgeOpacityDefault(50, style.name = style_name)
  setEdgeColorDefault("#AAAAAA", style.name = style_name)
  setEdgeLineWidthDefault(0.5, style.name = style_name)
  setBackgroundColorDefault("#FFFFFF", style.name = style_name)
  setVisualStyle(style_name)
  
  layoutNetwork("kamada-kawai")
  Sys.sleep(3)
  fitContent()
  
  exportImage(file.path(FIG_ALPACA, "FXS_ALPACA_network_final.png"), type = "PNG",
              resolution = 300, zoom = 200)
  exportNetwork(file.path(NETWORKS_DIR, "FXS_ALPACA_network.xgmml"), type = "xgmml")
  write.csv(edge_table, file.path(RESULTS_ALPACA, "FXS_ALPACA_edges.csv"), row.names = FALSE)
  write.csv(node_table, file.path(RESULTS_ALPACA, "FXS_ALPACA_nodes.csv"), row.names = FALSE)
  
  cat("Step 14 done. If labels remain hidden the network is dense —\n")
  cat("filter to fewer nodes or raise the font further.\n")
  
} else {
  # v5: Cytoscape-free export (works on Colab). Same module selection
  # as above; XGMML is written directly so Cytoscape can open it later.
  top_module_ids <- alpaca_df %>%
    group_by(Module) %>%
    summarise(mean_score = mean(Score), n = n(), .groups = "drop") %>%
    filter(n >= 5) %>% arrange(desc(mean_score)) %>% head(3) %>% pull(Module)
  top_nodes <- alpaca_df %>%
    filter(Module %in% top_module_ids) %>%
    dplyr::select(Gene, Module, Score, Type) %>%
    distinct(Gene, .keep_all = TRUE)
  edge_table <- edge_list_panda %>%
    filter(TF %in% top_nodes$Gene, Target %in% top_nodes$Gene) %>%
    dplyr::rename(Gene = Target)
  node_table <- top_nodes %>%
    dplyr::rename(id = Gene) %>%
    mutate(label = id, Module = as.character(Module), Score = round(Score, 4),
           NodeColor = ifelse(Type == "TF", "#E74C3C", "#3498DB"))
  cat("Top modules:", top_module_ids, "| nodes:", nrow(node_table),
      "| edges:", nrow(edge_table), "\n")

  xml_esc <- function(x) gsub("<", "&lt;", gsub(">", "&gt;",
                gsub("&", "&amp;", gsub("\"", "&quot;", as.character(x)))))
  set.seed(42)
  g_mod <- graph_from_data_frame(
    data.frame(from = edge_table$TF, to = edge_table$Gene, weight = edge_table$Score),
    directed = TRUE, vertices = node_table %>% dplyr::select(id, everything()))
  lay <- layout_with_kk(g_mod, weights = NA) * 120
  node_xml <- vapply(seq_len(vcount(g_mod)), function(i) {
    v <- V(g_mod)[i]
    sprintf(paste0('  <node id="%d" label="%s">\n',
                   '    <att name="Type" type="string" value="%s"/>\n',
                   '    <att name="Module" type="string" value="%s"/>\n',
                   '    <att name="Score" type="real" value="%s"/>\n',
                   '    <graphics type="%s" fill="%s" w="40" h="40" x="%.1f" y="%.1f"/>\n',
                   '  </node>'),
            i, xml_esc(v$name), xml_esc(v$Type), xml_esc(v$Module), v$Score,
            ifelse(v$Type == "TF", "ELLIPSE", "ROUND_RECTANGLE"), v$NodeColor,
            lay[i, 1], lay[i, 2])
  }, character(1))
  el <- as_edgelist(g_mod, names = FALSE)
  edge_xml <- sprintf(paste0('  <edge source="%d" target="%d" label="%s (pd) %s">\n',
                             '    <att name="weight" type="real" value="%s"/>\n',
                             '  </edge>'),
                      el[, 1], el[, 2], xml_esc(V(g_mod)$name[el[, 1]]),
                      xml_esc(V(g_mod)$name[el[, 2]]), E(g_mod)$weight)
  writeLines(c('<?xml version="1.0" encoding="UTF-8"?>',
               '<graph label="FXS_ALPACA_TopModules" directed="1"',
               '  xmlns="http://www.cs.rpi.edu/XGMML"',
               '  xmlns:xlink="http://www.w3.org/1999/xlink">',
               node_xml, edge_xml, '</graph>'),
             file.path(NETWORKS_DIR, "FXS_ALPACA_network.xgmml"))
  write_graph(g_mod, file.path(NETWORKS_DIR, "FXS_ALPACA_network.graphml"), format = "graphml")

  png(file.path(FIG_ALPACA, "FXS_ALPACA_network_final.png"),
      width = 3000, height = 3000, res = 300)
  plot(g_mod, layout = lay,
       vertex.color = V(g_mod)$NodeColor, vertex.frame.color = "white",
       vertex.shape = ifelse(V(g_mod)$Type == "TF", "circle", "square"),
       vertex.size = scales::rescale(V(g_mod)$Score, to = c(3, 9)),
       vertex.label = V(g_mod)$name, vertex.label.cex = 0.35,
       vertex.label.color = "black", edge.arrow.size = 0.1,
       edge.color = adjustcolor("grey60", alpha.f = 0.3), edge.width = 0.3,
       main = "FXS ALPACA top differential modules")
  dev.off()
  write.csv(edge_table, file.path(RESULTS_ALPACA, "FXS_ALPACA_edges.csv"), row.names = FALSE)
  write.csv(node_table, file.path(RESULTS_ALPACA, "FXS_ALPACA_nodes.csv"), row.names = FALSE)
  rm(g_mod, lay, node_xml, edge_xml, el)
  cat("Step 14 done (XGMML + GraphML + PNG written without Cytoscape)\n")
}

# ============================================================
# STEP 15 — VALIDATION COHORT: GSE218676
# ============================================================

cat("\n=== Step 15: Validation cohort — GSE218676 ===\n")

v_raw <- read.delim(gzfile(file.path(DATA_DIR, "GSE218676_raw_counts_GRCh38.p13_NCBI.tsv.gz")),
                    header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(colnames(v_raw)[1] == "GeneID")

v_counts <- as.matrix(v_raw[, -1])
rownames(v_counts) <- as.character(v_raw$GeneID)
storage.mode(v_counts) <- "integer"
v_gsm <- colnames(v_counts)
cat("GSE218676 — genes:", nrow(v_counts), "| samples:", ncol(v_counts), "\n")
cat("NOTE: NCBI reprocessed 20 of 26 submitted samples. FXS_426_NPC\n")
cat("(both replicates) is absent, giving 5 FXS vs 5 Control.\n")

v_meta_all <- read.csv(file.path(DATA_DIR, "GSE218676_metadata.csv"), row.names = 1,
                       stringsAsFactors = FALSE)
missing <- setdiff(v_gsm, rownames(v_meta_all))
if (length(missing) > 0)
  stop("GSMs absent from metadata: ", paste(missing, collapse = ", "))
v_meta <- v_meta_all[v_gsm, ]

v_meta$line      <- sub("_RNAseq_Rep[0-9]+$", "", v_meta$title)
v_meta$donor     <- sub("_(NPC|iPSC)$", "", v_meta$line)
v_meta$cell_type <- trimws(sub("^cell type:", "", v_meta$characteristics_ch1))
v_meta$genotype  <- case_when(
  grepl("^PMcut", v_meta$donor) ~ "Premutation_edited",
  grepl("^PM_",   v_meta$donor) ~ "Premutation",
  grepl("^FXS",   v_meta$donor) ~ "FXS",
  grepl("^NL",    v_meta$donor) ~ "Control",
  TRUE                          ~ "Unassigned")

cat("\n--- Validation cohort annotation ---\n")
print(data.frame(GSM = v_gsm, line = v_meta$line,
                 cell_type = v_meta$cell_type,
                 genotype = v_meta$genotype, row.names = NULL))
if (any(v_meta$genotype == "Unassigned")) stop("Unassigned genotypes")

v_keep <- rep(TRUE, nrow(v_meta))
if (VALID_EXCLUDE_PREMUTATION)
  v_keep <- v_keep & !grepl("^Premutation", v_meta$genotype)
if (VALID_MODE == "npc")
  v_keep <- v_keep & v_meta$cell_type == "iPSC derived NPC"

v_meta   <- v_meta[v_keep, ]
v_counts <- v_counts[, v_keep]
cat("\nSamples retained:", ncol(v_counts),
    "| lines:", length(unique(v_meta$line)),
    "| donors:", length(unique(v_meta$donor)), "\n")

v_design_tbl <- v_meta %>% distinct(line, donor, cell_type, genotype)
cat("\nCell type x genotype (lines):\n")
print(table(v_design_tbl$cell_type, v_design_tbl$genotype))

v_counts_f <- v_counts[rowSums(v_counts) >= 10, ]
cat("\nGenes after filtering:", nrow(v_counts_f), "\n")

v_coldata <- data.frame(
  condition = factor(v_meta$genotype, levels = c("Control", "FXS")),
  cell_type = factor(gsub("[^A-Za-z]", "_", v_meta$cell_type)),
  line      = factor(v_meta$line),
  donor     = factor(v_meta$donor),
  row.names = colnames(v_counts_f))

v_formula <- if (VALID_MODE == "pooled" && nlevels(v_coldata$cell_type) > 1)
  ~ cell_type + condition else ~ condition
cat("Validation design:", deparse(v_formula), "\n")

v_dds <- DESeqDataSetFromMatrix(v_counts_f, v_coldata, v_formula)
v_dds <- collapseReplicates(v_dds, groupby = v_dds$line)
colData(v_dds) <- droplevels(colData(v_dds))
cat("Biological units after collapsing:", ncol(v_dds), "\n")

v_post <- table(colData(v_dds)$cell_type, colData(v_dds)$condition)
cat("\nAfter collapsing:\n"); print(v_post)
if (any(v_post == 0))
  stop("Rank-deficient validation design — set VALID_MODE <- 'npc'")

v_dds <- DESeq(v_dds)
v_res <- results(v_dds, contrast = c("condition", "FXS", "Control"))
deg_valid <- as.data.frame(v_res)
deg_valid$GeneID <- rownames(deg_valid)
deg_valid$Symbol <- id2sym[deg_valid$GeneID]
deg_valid <- deg_valid[order(deg_valid$padj, na.last = TRUE), ]

cat("\nValidation FXS:", sum(colData(v_dds)$condition == "FXS"),
    "| Control:", sum(colData(v_dds)$condition == "Control"), "\n")
summary(v_res, alpha = 0.05)
write.csv(deg_valid, file.path(RESULTS_VALID, "GSE218676_validation_DEA_all.csv"))

sig_valid <- deg_valid %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(padj)
write.csv(sig_valid, file.path(RESULTS_VALID, "GSE218676_validation_DEA_significant.csv"))
cat("Validation DEGs:", nrow(sig_valid), "\n")

# v5: shrunken validation LFCs for the cross-cohort correlation
v_shr <- lfcShrink(v_dds, coef = "condition_FXS_vs_Control", type = "apeglm", quiet = TRUE)
deg_valid$log2FC_shrunk <- v_shr[deg_valid$GeneID, "log2FoldChange"]
rm(v_shr)

# ---- v5: donor-aware sensitivity analysis ----
# Two donors (FXS_421, NL_18) each contribute two lines (NPC + iPSC), so
# lines are not independent. DESeq2 has no random effects; limma-voom with
# duplicateCorrelation treats donor as a blocking factor and gives the
# effective-n-corrected answer. Agreement with DESeq2 is reported.
v_cd   <- as.data.frame(colData(v_dds))
v_dge  <- DGEList(counts(v_dds))
v_keep_g <- filterByExpr(v_dge, group = v_cd$condition)
v_dge  <- calcNormFactors(v_dge[v_keep_g, , keep.lib.sizes = FALSE])
v_des  <- model.matrix(if (nlevels(v_cd$cell_type) > 1) ~ cell_type + condition else ~ condition,
                       data = v_cd)
v_voom <- voom(v_dge, v_des)
v_corr <- duplicateCorrelation(v_voom, v_des, block = v_cd$donor)$consensus.correlation
v_voom <- voom(v_dge, v_des, block = v_cd$donor, correlation = v_corr)
v_corr <- duplicateCorrelation(v_voom, v_des, block = v_cd$donor)$consensus.correlation
v_lfit <- eBayes(lmFit(v_voom, v_des, block = v_cd$donor, correlation = v_corr))
v_tt   <- topTable(v_lfit, coef = "conditionFXS", number = Inf)
v_tt$GeneID <- rownames(v_tt); v_tt$Symbol <- id2sym[v_tt$GeneID]
write.csv(v_tt, file.path(RESULTS_VALID, "GSE218676_validation_limma_donorblock_all.csv"))

sig_valid_limma <- v_tt %>% filter(adj.P.Val < 0.05, abs(logFC) > 1)
both_tested <- intersect(v_tt$GeneID, deg_valid$GeneID[!is.na(deg_valid$log2FoldChange)])
donor_sens <- data.frame(
  Metric = c("Within-donor correlation (consensus)",
             "DEGs DESeq2 (lines as independent)",
             "DEGs limma-voom (donor blocked)",
             "DEGs shared by both", "LFC correlation DESeq2 vs limma (Pearson)",
             "FMR1 log2FC limma", "FMR1 FDR limma"),
  Value = c(round(v_corr, 3), nrow(sig_valid), nrow(sig_valid_limma),
            length(intersect(sig_valid$GeneID, sig_valid_limma$GeneID)),
            round(cor(v_tt[both_tested, "logFC"],
                      deg_valid[match(both_tested, deg_valid$GeneID), "log2FoldChange"]), 3),
            signif(v_tt$logFC[which(v_tt$Symbol == "FMR1")[1]], 4),
            signif(v_tt$adj.P.Val[which(v_tt$Symbol == "FMR1")[1]], 3)))
print(donor_sens)
write.csv(donor_sens, file.path(RESULTS_VALID, "validation_donor_sensitivity.csv"), row.names = FALSE)
rm(v_dge, v_voom, v_lfit, v_keep_g)

v_fmr1 <- deg_valid[which(deg_valid$Symbol == "FMR1"), ]
if (nrow(v_fmr1) > 0) {
  cat("\nFMR1 in validation cohort:\n")
  print(v_fmr1[, c("Symbol","baseMean","log2FoldChange","pvalue","padj")])
}

rm(v_raw, v_counts, v_counts_f); gc(verbose = FALSE)

# ---- 15D: OVERLAP ----

primary_genes <- unique(sig_degs$Symbol[!is.na(sig_degs$Symbol)])
valid_genes   <- unique(sig_valid$Symbol[!is.na(sig_valid$Symbol)])
cat("\nPrimary DEGs (GSE146339 organoids):", length(primary_genes), "\n")
cat("Validation DEGs (GSE218676):", length(valid_genes), "\n")

overlap_genes <- intersect(primary_genes, valid_genes)
cat("Overlapping DEGs:", length(overlap_genes), "\n")
print(overlap_genes)

if (length(overlap_genes) > 0) {
  p_lfc <- sapply(overlap_genes, function(g)
    mean(sig_degs$log2FoldChange[sig_degs$Symbol == g], na.rm = TRUE))
  v_lfc <- sapply(overlap_genes, function(g)
    mean(sig_valid$log2FoldChange[sig_valid$Symbol == g], na.rm = TRUE))
  
  dir_check2 <- data.frame(Gene = overlap_genes,
                           Primary_LFC = p_lfc, Validation_LFC = v_lfc) %>%
    filter(!is.na(Primary_LFC), !is.na(Validation_LFC))
  direction_concordant <- sum(sign(dir_check2$Primary_LFC) ==
                                sign(dir_check2$Validation_LFC))
  common_for_direction <- nrow(dir_check2)
  
  cat("\nDirection check:\n"); print(dir_check2)
  cat("Concordant:", direction_concordant, "/", common_for_direction, "\n")
  write.csv(dir_check2, file.path(RESULTS_VALID, "crossvalidation_overlap_genes.csv"), row.names = FALSE)
  
  dir_check2$Direction <- ifelse(sign(dir_check2$Primary_LFC) ==
                                   sign(dir_check2$Validation_LFC),
                                 "Concordant", "Discordant")
  ggplot(dir_check2, aes(x = Primary_LFC, y = Validation_LFC,
                         color = Direction, label = Gene)) +
    geom_point(size = 4, alpha = 0.9) +
    geom_text_repel(size = 3.5, fontface = "bold", max.overlaps = 15) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("Concordant" = "#2ECC71",
                                  "Discordant" = "#E74C3C"), name = "Direction") +
    labs(title = "Direction Concordance — Organoids vs NPC/B-cell",
         subtitle = paste0(direction_concordant, "/", common_for_direction,
                           " genes concordant"),
         x = "GSE146339 log2FC (organoids)",
         y = "GSE218676 log2FC (NPC/B-cell/iPSC)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  ggsave(file.path(FIG_VALID, "crossvalidation_direction_concordance.png"),
         width = 8, height = 6, dpi = 300)
} else {
  cat("No overlapping DEGs — the Step 16 rank-based analysis still applies\n")
  cat("and is the more informative test at small overlap.\n")
  dir_check2 <- data.frame(); direction_concordant <- 0; common_for_direction <- 0
}

overlap_in_panda  <- intersect(overlap_genes,
                               unique(c(edge_list_panda$TF, edge_list_panda$Target)))
overlap_in_alpaca <- intersect(overlap_genes, alpaca_df$Gene)
cat("Overlap in PANDA:", length(overlap_in_panda),
    "| in ALPACA:", length(overlap_in_alpaca), "\n")

ggvenn(list("Organoids\n(GSE146339)" = primary_genes,
            "NPC/B-cell\n(GSE218676)" = valid_genes),
       fill_color = c("#E74C3C", "#3498DB"), stroke_size = 0.5,
       set_name_size = 4) +
  labs(title = "DEG Overlap — Two Independent FXS Cohorts",
       subtitle = paste0(length(overlap_genes), " genes replicated")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11))
ggsave(file.path(FIG_VALID, "crossvalidation_venn.png"), width = 8, height = 6, dpi = 300)

cat("Step 15 done\n")

# ============================================================
# STEP 16 — CROSS-COHORT CONCORDANCE
# ============================================================

cat("\n=== STEP 16: Cross-cohort concordance ===\n")

primary_rank <- res_df %>%
  filter(!is.na(Symbol), Symbol != "", !is.na(pvalue), !is.na(log2FoldChange)) %>%
  mutate(rank_metric = -log10(pvalue) * sign(log2FoldChange)) %>%
  group_by(Symbol) %>% slice_max(abs(rank_metric), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(Gene = Symbol, log2FC = log2FoldChange, log2FC_shrunk, pvalue, rank_metric) %>%
  as.data.frame()
cat("Primary cohort testable genes:", nrow(primary_rank), "\n")

valid_rank <- deg_valid %>%
  filter(!is.na(Symbol), Symbol != "", !is.na(pvalue), !is.na(log2FoldChange)) %>%
  mutate(rank_metric = -log10(pvalue) * sign(log2FoldChange)) %>%
  group_by(Symbol) %>% slice_max(abs(rank_metric), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(Gene = Symbol, log2FC = log2FoldChange, log2FC_shrunk, pvalue, rank_metric) %>%
  as.data.frame()
cat("Validation cohort testable genes:", nrow(valid_rank), "\n")

shared_genes <- intersect(primary_rank$Gene, valid_rank$Gene)
cat("Genes measured in both cohorts (background):", length(shared_genes), "\n")

primary_shared <- primary_rank[match(shared_genes, primary_rank$Gene), ]
valid_shared   <- valid_rank[match(shared_genes, valid_rank$Gene), ]

# ---- 16B: OVERLAP SIGNIFICANCE ----

primary_sig_shared <- intersect(primary_genes, shared_genes)
valid_sig_shared   <- intersect(valid_genes,   shared_genes)
overlap_shared     <- intersect(primary_sig_shared, valid_sig_shared)

N <- length(shared_genes); K <- length(primary_sig_shared)
M <- length(valid_sig_shared); k <- length(overlap_shared)

cat("\n--- Overlap significance ---\n")
cat("Background N:", N, "| Primary DEGs:", K,
    "| Validation DEGs:", M, "| Overlap:", k, "\n")

expected_overlap <- if (N > 0) K * M / N else NA
fold_enrichment  <- if (!is.na(expected_overlap) && expected_overlap > 0)
  k / expected_overlap else NA
hyper_p <- phyper(k - 1, M, N - M, K, lower.tail = FALSE)
fisher_tab <- matrix(c(k, K - k, M - k, N - K - M + k), nrow = 2,
                     dimnames = list(c("Valid_DEG","Valid_notDEG"),
                                     c("Primary_DEG","Primary_notDEG")))
fisher_res <- fisher.test(fisher_tab, alternative = "greater")

cat("Expected by chance:", round(expected_overlap, 2),
    "| Fold enrichment:", round(fold_enrichment, 2), "\n")
cat("Hypergeometric p:", format.pval(hyper_p, digits = 3),
    "| Fisher p:", format.pval(fisher_res$p.value, digits = 3),
    "| OR:", round(fisher_res$estimate, 2), "\n")
print(fisher_tab)

write.csv(data.frame(Background_N = N, Primary_DEGs = K, Validation_DEGs = M,
                     Observed_overlap = k, Expected_overlap = expected_overlap,
                     Fold_enrichment = fold_enrichment,
                     Hypergeometric_p = hyper_p, Fisher_p = fisher_res$p.value,
                     Fisher_OR = as.numeric(fisher_res$estimate)),
          file.path(RESULTS_VALID, "crossvalidation_overlap_significance.csv"), row.names = FALSE)

# ---- 16C: RRHO ----

run_rrho_manual <- function(rank1, rank2, genes, step = NULL) {
  n <- length(genes)
  if (is.null(step)) step <- max(1, floor(n / 100))
  ord1 <- genes[order(rank1, decreasing = TRUE)]
  ord2 <- genes[order(rank2, decreasing = TRUE)]
  cuts <- seq(step, n - step, by = step)
  mat <- matrix(0, nrow = length(cuts), ncol = length(cuts),
                dimnames = list(cuts, cuts))
  for (a in seq_along(cuts)) {
    s1 <- ord1[1:cuts[a]]
    for (b in seq_along(cuts)) {
      s2 <- ord2[1:cuts[b]]
      ov <- length(intersect(s1, s2))
      p  <- phyper(ov - 1, length(s2), n - length(s2), length(s1),
                   lower.tail = FALSE)
      mat[a, b] <- -log10(max(p, .Machine$double.xmin))
    }
  }
  list(matrix = mat, cuts = cuts, n = n)
}

rrho2_ok <- requireNamespace("RRHO2", quietly = TRUE) &&
  exists("RRHO2_initialize", envir = asNamespace("RRHO2"), inherits = FALSE)
if (!rrho2_ok && requireNamespace("RRHO2", quietly = TRUE))
  cat("Installed RRHO2 lacks RRHO2_initialize — using the built-in RRHO instead.\n")
if (rrho2_ok) {
  cat("\nRunning RRHO2...\n")
  rrho_obj <- getFromNamespace("RRHO2_initialize", "RRHO2")(
    data.frame(Gene = primary_shared$Gene, Value = primary_shared$rank_metric),
    data.frame(Gene = valid_shared$Gene,   Value = valid_shared$rank_metric),
    labels = c("Organoids (GSE146339)", "NPC/B-cell (GSE218676)"),
    log10.ind = TRUE)
  png(file.path(FIG_VALID, "RRHO2_heatmap.png"), width = 2000, height = 2000, res = 250)
  getFromNamespace("RRHO2_heatmap", "RRHO2")(rrho_obj); dev.off()
  rrho_max <- max(rrho_obj$hypermat, na.rm = TRUE)
  cat("RRHO2 peak -log10(p):", round(rrho_max, 2), "\n")
} else {
  cat("\nRRHO2 not installed — using built-in manual RRHO.\n")
  rrho_manual <- run_rrho_manual(primary_shared$rank_metric,
                                 valid_shared$rank_metric, shared_genes)
  rrho_max <- max(rrho_manual$matrix, na.rm = TRUE)
  cat("Manual RRHO peak -log10(p):", round(rrho_max, 2), "\n")
  rrho_long <- reshape2::melt(rrho_manual$matrix)
  colnames(rrho_long) <- c("Primary_rank", "Validation_rank", "neglog10p")
  ggplot(rrho_long, aes(x = Primary_rank, y = Validation_rank, fill = neglog10p)) +
    geom_tile() +
    scale_fill_viridis_c(name = "-log10(p)", option = "inferno") +
    labs(title = "Rank-Rank Hypergeometric Overlap (RRHO)",
         subtitle = paste0("Peak -log10(p) = ", round(rrho_max, 2),
                           " | ", length(shared_genes), " jointly-measured genes"),
         x = "Organoid rank cutoff (GSE146339)",
         y = "NPC/B-cell rank cutoff (GSE218676)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  ggsave(file.path(FIG_VALID, "RRHO_manual_heatmap.png"), width = 8, height = 7, dpi = 300)
  write.csv(rrho_long, file.path(RESULTS_VALID, "RRHO_manual_matrix.csv"), row.names = FALSE)
}

# ---- 16D: GLOBAL EFFECT-SIZE CORRELATION ----

lfc_df <- data.frame(Gene = shared_genes,
                     Primary_LFC = primary_shared$log2FC,
                     Validation_LFC = valid_shared$log2FC,
                     Primary_LFC_shrunk = primary_shared$log2FC_shrunk,
                     Validation_LFC_shrunk = valid_shared$log2FC_shrunk) %>%
  filter(!is.na(Primary_LFC), !is.na(Validation_LFC))

# v5: the same statistics on apeglm-shrunken LFCs. Unshrunk LFCs of
# low-count genes are mostly noise and pull r toward 0 regardless of biology.
lfc_shr <- lfc_df %>% filter(!is.na(Primary_LFC_shrunk), !is.na(Validation_LFC_shrunk))
pearson_shr  <- cor.test(lfc_shr$Primary_LFC_shrunk, lfc_shr$Validation_LFC_shrunk)
spearman_shr <- cor.test(lfc_shr$Primary_LFC_shrunk, lfc_shr$Validation_LFC_shrunk,
                         method = "spearman", exact = FALSE)
lfc_shr_eff  <- lfc_shr %>% filter(abs(Primary_LFC_shrunk) > 0.1,
                                   abs(Validation_LFC_shrunk) > 0.1)
conc_shr_eff <- mean(sign(lfc_shr_eff$Primary_LFC_shrunk) ==
                       sign(lfc_shr_eff$Validation_LFC_shrunk))
cat("Shrunken-LFC Pearson r:", round(pearson_shr$estimate, 3),
    "| Spearman:", round(spearman_shr$estimate, 3),
    "| sign concordance where both |shrunk LFC|>0.1:",
    round(100 * conc_shr_eff, 1), "% of", nrow(lfc_shr_eff), "genes\n")

pearson_res  <- cor.test(lfc_df$Primary_LFC, lfc_df$Validation_LFC,
                         method = "pearson")
spearman_res <- cor.test(lfc_df$Primary_LFC, lfc_df$Validation_LFC,
                         method = "spearman", exact = FALSE)
global_concordance <- mean(sign(lfc_df$Primary_LFC) == sign(lfc_df$Validation_LFC))
lfc_strong <- lfc_df %>%
  filter(abs(Primary_LFC) > 0.5 | abs(Validation_LFC) > 0.5)

cat("\n--- Global effect-size correlation ---\n")
cat("Genes compared:", nrow(lfc_df), "\n")
cat("Pearson r:", round(pearson_res$estimate, 3),
    "p:", format.pval(pearson_res$p.value, digits = 3), "\n")
cat("Spearman rho:", round(spearman_res$estimate, 3),
    "p:", format.pval(spearman_res$p.value, digits = 3), "\n")
cat("Genome-wide sign concordance:", round(global_concordance * 100, 1), "%\n")
cat("INTERPRET: concordance near 50% and r near 0 means the cohorts do\n")
cat("NOT agree transcriptome-wide. A significant DEG overlap can still\n")
cat("be enriched above chance. Report both — different claims.\n")
cat("Tissue differs between cohorts, so partial disagreement is expected.\n")

lfc_df$Highlight <- lfc_df$Gene %in% overlap_shared
ggplot(lfc_df, aes(x = Primary_LFC, y = Validation_LFC)) +
  geom_point(alpha = 0.15, size = 0.7, color = "grey50") +
  geom_point(data = subset(lfc_df, Highlight), color = "#E74C3C", size = 2.5) +
  geom_text_repel(data = subset(lfc_df, Highlight), aes(label = Gene),
                  size = 3, fontface = "bold", max.overlaps = 20) +
  geom_smooth(method = "lm", se = TRUE, color = "#2C3E50", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  labs(title = "Transcriptome-wide Effect-Size Concordance",
       subtitle = paste0("Pearson r = ", round(pearson_res$estimate, 3),
                         ", p = ", format.pval(pearson_res$p.value, digits = 3),
                         " | n = ", nrow(lfc_df), " shared genes"),
       x = "GSE146339 log2FC (organoids)",
       y = "GSE218676 log2FC (NPC/B-cell/iPSC)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
ggsave(file.path(FIG_VALID, "crossvalidation_global_LFC_correlation.png"), width = 8, height = 7, dpi = 300)
write.csv(lfc_df, file.path(RESULTS_VALID, "crossvalidation_shared_gene_LFC.csv"), row.names = FALSE)

write.csv(data.frame(
  Metric = c("Shared genes","Pearson r","Pearson p","Spearman rho","Spearman p",
             "Genome-wide sign concordance (%)","Sign concordance |LFC|>0.5 (%)",
             "RRHO peak -log10(p)","Overlap fold enrichment",
             "Overlap hypergeometric p",
             "Pearson r (shrunken LFC)","Pearson p (shrunken LFC)",
             "Spearman rho (shrunken LFC)",
             "Sign concordance, both |shrunk LFC|>0.1 (%)",
             "Genes with both |shrunk LFC|>0.1"),
  Value = c(nrow(lfc_df), round(pearson_res$estimate, 4), pearson_res$p.value,
            round(spearman_res$estimate, 4), spearman_res$p.value,
            round(global_concordance * 100, 1),
            round(mean(sign(lfc_strong$Primary_LFC) ==
                         sign(lfc_strong$Validation_LFC)) * 100, 1),
            round(rrho_max, 3), round(fold_enrichment, 3), hyper_p,
            round(pearson_shr$estimate, 4), pearson_shr$p.value,
            round(spearman_shr$estimate, 4),
            round(100 * conc_shr_eff, 1), nrow(lfc_shr_eff))),
  file.path(RESULTS_VALID, "crossvalidation_concordance_summary.csv"), row.names = FALSE)

cat("Step 16 done\n")

# ============================================================
# CHECKPOINT
# ============================================================

rm(v_dds, v_res, module_go_results, top_go_strict, top_go_relaxed,
   top_kegg_strict, top_kegg_relaxed, best_go, best_kegg,
   g_genie3, g_panda, genie3_filtered, panda_filtered, lioness_df)
gc(verbose = FALSE)

save.image(file.path(CHECKPOINTS_DIR, "checkpoint_step16.RData"))
cat("\n=== PART 2 COMPLETE — checkpoint_step16.RData saved ===\n")
cat("Restart R, then run part3_steps17-19.R\n")