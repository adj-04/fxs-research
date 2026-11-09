# ============================================================
# FXS PIPELINE v5 — STEP 17B: PANDA randomized-prior null
# Run as its own Rscript call after part 2. Loads checkpoint_step16.
#
# v5 changes
#  * Pooled null. v4 gave each edge an empirical p of
#    (#runs with null >= obs + 1) / 21, so the smallest possible p was
#    0.048 and no edge could ever pass BH — "0 at FDR" was guaranteed
#    by the design, not evidence of absent signal. v5 compares each
#    observed edge with the pooled distribution of every null edge
#    from every run (resolution ~ 1 / (runs x edges)).
#  * Resumable: state is saved after every run, so a Colab disconnect
#    costs at most one PANDA run. Rerun the same command to continue.
# ============================================================

library(pandaR)
library(igraph)
library(dplyr)
library(ggplot2)

source(file.path(Sys.getenv("FXS_REPO", "."), "scripts", "00_config.R"))
load(ckpt("checkpoint_step16.RData"))
source(file.path(Sys.getenv("FXS_REPO", "."), "scripts", "00_config.R"))
log_session("part17b")
cat("Checkpoint loaded.\n")

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

obs_edge   <- paste(edge_list_panda$TF, edge_list_panda$Target, sep = "__")
obs_scores <- edge_list_panda$Score

state_file <- ckpt("panda_null_state.rds")
if (file.exists(state_file)) {
  st <- readRDS(state_file)
  cat("Resuming PANDA null from run", st$done + 1, "/", N_PERM_PANDA, "\n")
} else {
  st <- list(done = 0,
             hist = make_null_hist(lo = -60, hi = 60, log = FALSE),
             exceed = integer(length(obs_scores)),   # v4-style per-edge counts
             sample = numeric(0), n_ok = 0)
}

t0 <- Sys.time()
while (st$done < N_PERM_PANDA) {
  p <- st$done + 1
  set.seed(1000 + p)                     # reproducible per run, resume-safe
  cat("  PANDA null run", p, "/", N_PERM_PANDA,
      "| elapsed", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

  motif_rand <- randomize_motif_degpres(motif_edgelist)
  motif_rand <- motif_rand[motif_rand$Gene %in% rownames(expr_final) &
                             motif_rand$TF %in% rownames(expr_final), ]
  pr <- if (nrow(motif_rand) == 0) NULL else tryCatch(
    pandaR::panda(motif = motif_rand, expr = as.data.frame(expr_final),
                  ppi = ppi_edgelist, alpha = 0.1, hamming = 0.001,
                  iter = 1000, output = "regulatory", zScale = TRUE,
                  progress = FALSE, randomize = "None",
                  cor.method = "pearson", scale.by.present = FALSE,
                  edgelist = TRUE, remove.missing.ppi = FALSE,
                  remove.missing.motif = FALSE, remove.missing.genes = FALSE,
                  mode = "union"),
    error = function(e) { cat("    failed:", conditionMessage(e), "\n"); NULL })

  if (!is.null(pr)) {
    rn <- slot(pr, "regNet")
    null_scores <- as.numeric(as.character(rn$Score))
    null_ids    <- paste(as.character(rn$TF), as.character(rn$Gene), sep = "__")
    st$hist   <- add_to_null_hist(st$hist, null_scores)
    st$sample <- c(st$sample, sample(null_scores, min(50000, length(null_scores))))
    m <- match(obs_edge, null_ids)
    hit <- !is.na(m)
    st$exceed[hit] <- st$exceed[hit] + (null_scores[m[hit]] >= obs_scores[hit])
    st$n_ok <- st$n_ok + 1
    rm(pr, rn, null_scores, null_ids, m, hit)
  }
  st$done <- p
  saveRDS(st, state_file)
  gc(verbose = FALSE)
}
cat("PANDA null complete:", st$n_ok, "successful runs\n")

panda_null_q95 <- null_hist_quantile(st$hist, 0.95)
edge_list_panda_sig <- edge_list_panda %>%
  mutate(Edge          = obs_edge,
         emp_p_pooled  = null_hist_p(st$hist, Score),
         emp_FDR_pooled = p.adjust(emp_p_pooled, method = "BH"),
         emp_p_peredge = (st$exceed + 1) / (st$n_ok + 1),       # v4 statistic
         emp_FDR_peredge = p.adjust(emp_p_peredge, method = "BH"),
         above_null95  = Score >= panda_null_q95) %>%
  arrange(emp_p_pooled, desc(Score))

panda_null_summary <- data.frame(
  Metric = c("Null runs", "Null edges pooled", "Null q95", "Observed 90th pct",
             "Smallest attainable p (per-edge, v4)", "Smallest attainable p (pooled, v5)",
             "Edges pooled p < 0.05", "Edges pooled FDR < 0.05",
             "Edges per-edge FDR < 0.05", "Edges above null q95"),
  Value = c(st$n_ok, st$hist$n, signif(panda_null_q95, 4), signif(panda_threshold, 4),
            signif(1 / (st$n_ok + 1), 3), signif(1 / (st$hist$n + 1), 3),
            sum(edge_list_panda_sig$emp_p_pooled < 0.05),
            sum(edge_list_panda_sig$emp_FDR_pooled < 0.05),
            sum(edge_list_panda_sig$emp_FDR_peredge < 0.05),
            sum(edge_list_panda_sig$above_null95)))
print(panda_null_summary)
write.csv(panda_null_summary, file.path(RESULTS_PANDA, "PANDA_null_summary.csv"), row.names = FALSE)
write.csv(edge_list_panda_sig, file.path(RESULTS_PANDA, "PANDA_edge_list_with_null.csv"), row.names = FALSE)

panda_null_pool <- st$sample
save(edge_list_panda_sig, panda_null_pool, panda_null_q95, panda_null_summary,
     file = ckpt("checkpoint_17B_only.RData"))

panda_cmp <- rbind(
  data.frame(Score = sample(panda_null_pool, min(100000, length(panda_null_pool))),
             Source = "Randomized prior"),
  data.frame(Score = sample(edge_list_panda$Score, min(100000, nrow(edge_list_panda))),
             Source = "Observed"))
ggplot(panda_cmp, aes(x = Score, fill = Source)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = panda_null_q95, linetype = "dashed") +
  scale_fill_manual(values = c("Observed" = "#E74C3C", "Randomized prior" = "#95A5A6")) +
  labs(title = "PANDA Edge Scores — Observed vs Randomized Prior",
       subtitle = paste0(st$n_ok, " degree-preserving randomizations | pooled null"),
       x = "PANDA edge score (z-scaled)", y = "Density") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
ggsave(file.path(FIG_PANDA, "PANDA_null_distribution.png"), width = 9, height = 6, dpi = 300)

cat("\n=== STEP 17B COMPLETE ===\n")
