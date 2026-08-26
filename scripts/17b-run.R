# ============================================================
# STEP 17B — PANDA RANDOMIZED-PRIOR NULL (standalone run)
# Run in a FRESH R session. Loads checkpoint_step16.RData only —
# lighter than step19, avoids carrying GENIE3/ALPACA null objects
# in memory while PANDA reruns 20x.
# ============================================================

library(pandaR)
library(igraph)
library(dplyr)
library(ggplot2)

setwd("~/Desktop/final report")
load("checkpoint_step16.RData")
cat("Checkpoint loaded. Objects:", length(ls()), "\n")

N_PERM_PANDA   <- 20
N_CORES        <- 2
set.seed(42)

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
                       paste(edge_list_panda$TF, edge_list_panda$Target, sep = "__"))
panda_exceed    <- setNames(rep(0L, length(obs_scores)), names(obs_scores))
panda_null_pool <- numeric(0)

t0 <- Sys.time()
for (p in seq_len(N_PERM_PANDA)) {
  cat("  PANDA null run", p, "/", N_PERM_PANDA,
      "| elapsed", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
  
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
                         paste(as.character(rn$TF), as.character(rn$Gene), sep = "__"))
  panda_null_pool <- c(panda_null_pool,
                       sample(null_edges, min(50000, length(null_edges))))
  common <- intersect(names(obs_scores), names(null_edges))
  panda_exceed[common] <- panda_exceed[common] +
    (null_edges[common] >= obs_scores[common])
  
  rm(pr, rn, null_edges); gc(verbose = FALSE)   # release before next run
}
cat("PANDA null complete in",
    round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

# ---- Extract results and save immediately (before plotting) ----

panda_emp_p    <- (panda_exceed + 1) / (N_PERM_PANDA + 1)
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
cat("Edges with empirical FDR < 0.05:",
    sum(edge_list_panda_sig$emp_FDR < 0.05, na.rm = TRUE), "\n")
cat("INTERPRET: if the null q95 sits close to the observed 90th\n")
cat("percentile, PANDA is largely reproducing its prior rather than\n")
cat("adding expression-driven signal. Report this in Limitations.\n")

write.csv(edge_list_panda_sig, "PANDA_edge_list_with_null.csv", row.names = FALSE)

save(edge_list_panda_sig, panda_null_pool, panda_null_q95, panda_emp_p,
     file = "checkpoint_17B_only.RData")
cat("Saved checkpoint_17B_only.RData\n")

# ---- Plot ----

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

cat("\n=== STEP 17B COMPLETE ===\n")