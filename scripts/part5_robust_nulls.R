# ============================================================
# FXS PIPELINE v5.1 — PART 5 : ROBUST NULL MODELS
#
# Replaces three tests from Part 3 / 17B that a review showed to be
# invalid or uninterpretable:
#
#   GENIE3  The v5 null used 250 trees vs 1,000 for the observed network,
#           which miscalibrates it (only 2% of edges at p < 0.05). Here the
#           null uses the same number of trees as the observed run, and a
#           leave-one-out check measures calibration directly (the fraction
#           of null edges with p < 0.05 should be ~5%).
#
#   PANDA   (a) Like-for-like randomized-prior comparison: ALL observed
#           edges vs ALL null edges (v5 compared positive observed edges
#           with all null edges). (b) New shuffled-expression null: real
#           motif + PPI prior with each gene's expression permuted across
#           samples. This directly asks whether the expression data move
#           PANDA edges beyond what random expression would.
#
#   ALPACA  ALPACA normalises node scores to sum to 2 within each module,
#           so module means are 2/size and cannot be compared with a pooled
#           node null. ALPACA was also run with FXS as the reference
#           network. Here: control is the reference, modules are scored by
#           their raw differential modularity Q_m, and significance comes
#           from FXS/control label permutations within timepoint (400
#           possible relabelings), with a max-statistic (family-wise)
#           correction so module selection is accounted for.
#
#   Local centrality of known TFs is recomputed within each corrected
#   module, as a percentile among TFs only.
#
# Requires: checkpoint_step16.RData (Part 2) and, for PANDA (a),
#           checkpoint panda_null_state.rds (17b-run.R).
# Every loop is resumable: rerun the same command after a disconnect.
#
# Stages (comma-separated, default all):
#   FXS_PART5_STAGES="genie3,panda,alpaca,summary"
# Sizes:
#   FXS_P5_GENIE3_PERM   (5)     FXS_P5_GENIE3_TREES (= observed, 1000)
#   FXS_P5_PANDA_SHUF    (5)
#   FXS_P5_ALPACA_PERM   (100, max 399; extend later, it resumes)
# ============================================================

suppressPackageStartupMessages({
  library(GENIE3); library(pandaR); library(netZooR); library(igraph)
  library(dplyr); library(ggplot2); library(data.table)
})

source(file.path(Sys.getenv("FXS_REPO", "."), "scripts", "00_config.R"))
load(ckpt("checkpoint_step16.RData"))
source(file.path(Sys.getenv("FXS_REPO", "."), "scripts", "00_config.R"))
log_session("part5")
cat("Checkpoint loaded.\n")

STAGES <- trimws(strsplit(Sys.getenv("FXS_PART5_STAGES",
                                     "genie3,panda,alpaca,summary"), ",")[[1]])
OUT  <- file.path("results", "robust_nulls")
FIGO <- file.path("figures", "robust_nulls")
dir.create(OUT,  recursive = TRUE, showWarnings = FALSE)
dir.create(FIGO, recursive = TRUE, showWarnings = FALSE)
cat("Stages:", paste(STAGES, collapse = ", "), "\n")

ts <- function() format(Sys.time(), "%H:%M:%S")
hist_with_counts <- function(template, counts, n) { template$counts <- counts; template$n <- n; template }

# ============================================================
# STAGE 1 — GENIE3: matched-tree null with calibration check
# ============================================================
if ("genie3" %in% STAGES) {
  cat("\n=== GENIE3 null (matched trees) ===\n")
  P_G   <- env_num("FXS_P5_GENIE3_PERM", 5)
  TREES <- env_num("FXS_P5_GENIE3_TREES", N_TREES_OBS)
  cat("Permutations:", P_G, "| trees per target:", TREES,
      "(observed network used", N_TREES_OBS, ")\n")

  stopifnot(all(rownames(genie3_weights) %in% tf_in_data))
  self_mask <- outer(rownames(genie3_weights), colnames(genie3_weights), "==")
  h0 <- make_null_hist(lo = -12, hi = 0, log = TRUE)

  sf <- ckpt("p5_genie3_state.rds")
  st <- if (file.exists(sf)) readRDS(sf) else NULL
  if (!is.null(st) && st$trees != TREES) {
    cat("Saved state used", st$trees, "trees — starting fresh.\n"); st <- NULL
  }
  if (is.null(st)) st <- list(done = 0, trees = TREES, counts = list(), n = list(), samples = list())
  if (st$done > 0) cat("Resuming after permutation", st$done, "\n")

  permute_rows <- function(m) { o <- t(apply(m, 1, sample)); dimnames(o) <- dimnames(m); o }
  t0 <- Sys.time(); start_done <- st$done
  while (st$done < P_G) {
    p <- st$done + 1
    set.seed(5000 + p)
    w <- GENIE3(permute_rows(expr_subset), regulators = tf_in_data,
                nTrees = TREES, nCores = N_CORES, verbose = FALSE)
    w <- w[rownames(genie3_weights), colnames(genie3_weights), drop = FALSE]
    wv <- w[!self_mask]; wv <- wv[wv > 0]
    hh <- add_to_null_hist(h0, wv)
    st$counts[[p]]  <- hh$counts
    st$n[[p]]       <- hh$n
    st$samples[[p]] <- sample(wv, min(2e5, length(wv)))
    st$done <- p
    saveRDS(st, sf)
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("[%s] GENIE3 permutation %d/%d | %.1f min | ETA %.0f min\n", ts(), p, P_G,
                el, el / (p - start_done) * (P_G - p)))
    rm(w, wv, hh); gc(verbose = FALSE)
  }

  h_all <- hist_with_counts(h0, Reduce(`+`, st$counts), sum(unlist(st$n)))
  obs_p <- null_hist_p(h_all, edge_list_genie3$Weight)
  obs_fdr <- p.adjust(obs_p, method = "BH")

  # Leave-one-out calibration: p-values of held-out null edges
  calib <- do.call(rbind, lapply(seq_len(st$done), function(k) {
    hk <- hist_with_counts(h0, h_all$counts - st$counts[[k]], h_all$n - st$n[[k]])
    pk <- null_hist_p(hk, st$samples[[k]])
    data.frame(heldout = k, frac_p05 = mean(pk < 0.05), frac_p01 = mean(pk < 0.01))
  }))
  cat("Calibration (held-out null edges; expect ~0.05 and ~0.01):\n"); print(calib)

  g3_sum <- data.frame(
    Metric = c("Permutations", "Trees per target (null)", "Trees per target (observed)",
               "Edges tested", "Null values pooled", "Smallest attainable p",
               "Calibration: null fraction p<0.05 (mean)", "Calibration: null fraction p<0.01 (mean)",
               "Observed fraction p<0.05", "Observed edges FDR<0.05", "Observed edges FDR<0.10"),
    Value = c(st$done, TREES, N_TREES_OBS, length(obs_p), h_all$n, signif(1 / (h_all$n + 1), 3),
              round(mean(calib$frac_p05), 4), round(mean(calib$frac_p01), 4),
              round(mean(obs_p < 0.05), 4), sum(obs_fdr < 0.05), sum(obs_fdr < 0.10)))
  print(g3_sum, row.names = FALSE)
  write.csv(g3_sum, file.path(OUT, "GENIE3_matched_null_summary.csv"), row.names = FALSE)
  write.csv(calib,  file.path(OUT, "GENIE3_matched_null_calibration.csv"), row.names = FALSE)

  g3_edges <- data.frame(edge_list_genie3[, c("TF", "Target", "Weight")], p = obs_p, FDR = obs_fdr)
  keep <- g3_edges$FDR < 0.10 | rank(g3_edges$p, ties.method = "first") <= 10000
  fwrite(g3_edges[keep, ][order(g3_edges$p[keep]), ],
         file.path(OUT, "GENIE3_matched_null_top_edges.csv"))
  tf_sig <- g3_edges %>% filter(FDR < 0.05) %>% count(TF, name = "n_FDR05_edges") %>%
    arrange(desc(n_FDR05_edges))
  write.csv(tf_sig, file.path(OUT, "GENIE3_matched_null_edges_per_TF.csv"), row.names = FALSE)
  cat("Known TFs — edges at FDR<0.05:\n")
  print(tf_sig[tf_sig$TF %in% known_tf, ])

  png(file.path(FIGO, "GENIE3_matched_null_pvalues.png"), width = 2000, height = 900, res = 200)
  par(mfrow = c(1, 2))
  hist(sample(obs_p, min(1e6, length(obs_p))), breaks = 50, col = "#3498DB", border = "white",
       main = "Observed GENIE3 edges", xlab = "Pooled-null p-value")
  hk <- hist_with_counts(h0, h_all$counts - st$counts[[1]], h_all$n - st$n[[1]])
  hist(null_hist_p(hk, st$samples[[1]]), breaks = 50, col = "#95A5A6", border = "white",
       main = "Held-out null edges (should be flat)", xlab = "Pooled-null p-value")
  dev.off()
  rm(g3_edges, obs_p, obs_fdr, self_mask); gc(verbose = FALSE)
}

# ============================================================
# STAGE 2 — PANDA: like-for-like prior null + shuffled-expression null
# ============================================================
if ("panda" %in% STAGES) {
  cat("\n=== PANDA ===\n")
  run_panda <- function(motif, expr) {
    pandaR::panda(motif = motif, expr = as.data.frame(expr), ppi = ppi_edgelist,
                  alpha = 0.1, hamming = 0.001, iter = 1000,
                  output = "regulatory", zScale = TRUE, progress = FALSE,
                  randomize = "None", cor.method = "pearson", scale.by.present = FALSE,
                  edgelist = TRUE, remove.missing.ppi = FALSE,
                  remove.missing.motif = FALSE, remove.missing.genes = FALSE, mode = "union")
  }
  regnet <- function(pr) {
    rn <- slot(pr, "regNet")
    list(id = paste(as.character(rn$TF), as.character(rn$Gene), sep = "__"),
         score = as.numeric(as.character(rn$Score)))
  }

  # ---- observed network, ALL edges (rerun; deterministic) ----
  obs_f <- ckpt("p5_panda_obs.rds")
  if (file.exists(obs_f)) { obs <- readRDS(obs_f); cat("Observed PANDA loaded.\n") } else {
    cat("[", ts(), "] Rerunning observed PANDA to keep all edges...\n")
    obs <- regnet(run_panda(motif_edgelist, expr_final)); saveRDS(obs, obs_f)
  }
  chk_id <- paste(edge_list_panda$TF, edge_list_panda$Target, sep = "__")
  chk <- cor(edge_list_panda$Score, obs$score[match(chk_id, obs$id)], use = "complete.obs")
  cat("Observed edges:", length(obs$score), "| agreement with Part 1 (positive edges), r =",
      round(chk, 6), "\n")

  # ---- (a) randomized-prior null, like for like ----
  prior_rows <- NULL
  pn_f <- ckpt("panda_null_state.rds")
  if (file.exists(pn_f)) {
    pn <- readRDS(pn_f); hP <- pn$hist
    p_obs <- null_hist_p(hP, obs$score)
    q95 <- null_hist_quantile(hP, 0.95); q05 <- null_hist_quantile(hP, 0.05)
    grid <- seq(-8, 15, by = 0.01)
    cdf_null <- cumsum(hP$counts)[.hist_bin(hP, grid)] / hP$n
    cdf_obs  <- ecdf(obs$score)(grid)
    ks_D <- max(abs(cdf_obs - cdf_null))
    qs <- c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
    qtab <- data.frame(quantile = qs,
                       observed = quantile(obs$score, qs),
                       null = vapply(qs, function(q) null_hist_quantile(hP, q), numeric(1)))
    print(qtab, row.names = FALSE)
    prior_rows <- data.frame(
      Metric = c("Prior-null runs", "Observed edges (all)", "Null q95 (all edges)",
                 "Observed fraction above null q95", "Observed fraction below null q05",
                 "Max CDF difference (KS D)", "Observed edges pooled FDR<0.05 (upper tail)"),
      Value = c(pn$n_ok, length(obs$score), signif(q95, 4),
                round(mean(obs$score > q95), 4), round(mean(obs$score < q05), 4),
                round(ks_D, 4), sum(p.adjust(p_obs, "BH") < 0.05)))
    print(prior_rows, row.names = FALSE)
    write.csv(qtab, file.path(OUT, "PANDA_prior_null_quantiles.csv"), row.names = FALSE)

    set.seed(11)
    dd <- rbind(data.frame(Score = sample(obs$score, 2e5), Source = "Observed (all edges)"),
                data.frame(Score = sample(pn$sample, min(2e5, length(pn$sample))),
                           Source = "Randomized prior (all edges)"))
    ggplot(dd, aes(Score, fill = Source)) + geom_density(alpha = 0.5) +
      geom_vline(xintercept = q95, linetype = "dashed") +
      scale_fill_manual(values = c("#E74C3C", "#95A5A6")) +
      labs(title = "PANDA edge scores: observed vs degree-preserving randomized prior",
           subtitle = sprintf("All edges in both | %d randomizations | KS D = %.3f", pn$n_ok, ks_D),
           x = "PANDA edge score (z-scaled)", y = "Density") + theme_bw()
    ggsave(file.path(FIGO, "PANDA_prior_null_all_edges.png"), width = 9, height = 5.5, dpi = 300)
    rm(pn, dd, p_obs); gc(verbose = FALSE)
  } else cat("panda_null_state.rds not found — skipping prior-null comparison.\n")

  # ---- (b) shuffled-expression null ----
  K <- env_num("FXS_P5_PANDA_SHUF", 5)
  for (k in seq_len(K)) {
    f <- ckpt(sprintf("p5_panda_shuf_%02d.rds", k))
    if (file.exists(f)) next
    cat(sprintf("[%s] PANDA with shuffled expression %d/%d\n", ts(), k, K))
    set.seed(6000 + k)
    ex <- t(apply(expr_final, 1, sample)); dimnames(ex) <- dimnames(expr_final)
    r <- regnet(run_panda(motif_edgelist, ex))
    saveRDS(r$score[match(obs$id, r$id)], f)
    rm(ex, r); gc(verbose = FALSE)
  }
  S <- sapply(seq_len(K), function(k) readRDS(ckpt(sprintf("p5_panda_shuf_%02d.rds", k))))
  ok <- stats::complete.cases(S) & is.finite(obs$score)
  S <- S[ok, , drop = FALSE]; o <- obs$score[ok]; ids <- obs$id[ok]

  prior_vec <- numeric(length(ids))
  mp <- match(ids, paste(motif_edgelist$TF, motif_edgelist$Gene, sep = "__"))
  prior_vec[!is.na(mp)] <- motif_edgelist$Score[mp[!is.na(mp)]]

  r_os <- apply(S, 2, function(s) cor(o, s))
  r_ss <- combn(K, 2, function(ix) cor(S[, ix[1]], S[, ix[2]]))
  r_prior_obs  <- cor(o, prior_vec)
  r_prior_shuf <- mean(apply(S, 2, function(s) cor(s, prior_vec)))

  # Edge-level: does real expression move an edge beyond shuffled-expression spread?
  scale_k <- sqrt((1 + 1 / K) / (1 + 1 / (K - 1)))
  d_obs <- o - rowMeans(S)
  hD0 <- make_null_hist(lo = 0, hi = 40, log = FALSE)
  d0_counts <- vector("list", K); d0_n <- numeric(K); d0_samp <- vector("list", K)
  for (k in seq_len(K)) {
    d0 <- abs((S[, k] - rowMeans(S[, -k, drop = FALSE])) * scale_k)
    hk <- add_to_null_hist(hD0, d0)
    d0_counts[[k]] <- hk$counts; d0_n[k] <- hk$n
    d0_samp[[k]] <- sample(d0, min(2e5, length(d0)))
  }
  hD <- hist_with_counts(hD0, Reduce(`+`, d0_counts), sum(d0_n))
  p_d <- null_hist_p(hD, abs(d_obs)); fdr_d <- p.adjust(p_d, "BH")
  calib_d <- sapply(seq_len(K), function(k) {
    hk <- hist_with_counts(hD0, hD$counts - d0_counts[[k]], hD$n - d0_n[k])
    mean(null_hist_p(hk, d0_samp[[k]]) < 0.05)
  })

  expr_rows <- data.frame(
    Metric = c("Shuffled-expression runs",
               "Correlation observed vs shuffled-expression network (mean)",
               "Correlation between shuffled-expression networks (mean)",
               "Correlation with motif prior: observed",
               "Correlation with motif prior: shuffled expression (mean)",
               "Mean |observed - shuffled mean|",
               "Calibration: held-out fraction p<0.05 (mean)",
               "Edges moved by real expression, FDR<0.05",
               "  of which higher with real expression",
               "Edges moved by real expression, FDR<0.10"),
    Value = c(K, round(mean(r_os), 4), round(mean(r_ss), 4),
              round(r_prior_obs, 4), round(r_prior_shuf, 4),
              round(mean(abs(d_obs)), 4), round(mean(calib_d), 4),
              sum(fdr_d < 0.05), sum(fdr_d < 0.05 & d_obs > 0), sum(fdr_d < 0.10)))
  print(expr_rows, row.names = FALSE)

  pe <- data.frame(Edge = ids, TF = sub("__.*", "", ids), Target = sub(".*__", "", ids),
                   Observed = o, Shuffled_mean = rowMeans(S), Diff = d_obs, p = p_d, FDR = fdr_d)
  pe <- pe[order(pe$p), ]
  fwrite(head(pe[pe$FDR < 0.10 | seq_len(nrow(pe)) <= 10000, ], 50000),
         file.path(OUT, "PANDA_expression_effect_top_edges.csv"))
  tf_expr <- pe %>% filter(FDR < 0.05) %>% count(TF, name = "n_FDR05_edges") %>% arrange(desc(n_FDR05_edges))
  write.csv(tf_expr, file.path(OUT, "PANDA_expression_effect_edges_per_TF.csv"), row.names = FALSE)
  cat("Known TFs — edges moved by expression at FDR<0.05:\n"); print(tf_expr[tf_expr$TF %in% known_tf, ])

  write.csv(rbind(prior_rows, expr_rows), file.path(OUT, "PANDA_robust_null_summary.csv"), row.names = FALSE)

  set.seed(12); ix <- sample(length(o), min(2e5, length(o)))
  png(file.path(FIGO, "PANDA_shuffled_expression.png"), width = 2000, height = 900, res = 200)
  par(mfrow = c(1, 2))
  plot(rowMeans(S)[ix], o[ix], pch = ".", col = adjustcolor("#2C3E50", 0.3),
       xlab = "Shuffled-expression network (mean of runs)", ylab = "Observed network",
       main = sprintf("PANDA edges (r = %.3f; shuffle-vs-shuffle r = %.3f)", mean(r_os), mean(r_ss)))
  abline(0, 1, col = "#E74C3C")
  hist(p_d, breaks = 50, col = "#3498DB", border = "white",
       main = "Edge-level expression effect", xlab = "p-value (|obs - shuffled| vs held-out shuffles)")
  dev.off()
  rm(S, o, ids, pe, d_obs, p_d, fdr_d, obs, prior_vec); gc(verbose = FALSE)
}

# ============================================================
# STAGE 3 — ALPACA: control as reference, label-permutation null
# ============================================================
if ("alpaca" %in% STAGES) {
  cat("\n=== ALPACA (control reference, label permutation) ===\n")
  e_tf   <- sub("__.*", "", rownames(lioness_mat))
  e_gene <- sub(".*__", "", rownames(lioness_mat))
  samp <- colnames(lioness_mat)
  tp   <- as.character(timepoint[samp])
  cond <- as.character(condition[samp])
  stopifnot(!anyNA(tp), !anyNA(cond))
  print(table(timepoint = tp, condition = cond))

  # Same logic as netZooR::alpaca(), but returns raw module contributions.
  # (alpaca() normalises node scores within each module, so module totals
  # are always 2 and carry no information.)
  run_alpaca <- function(ctrl_w, cond_w, seed, keep_dwbm = FALSE) {
    set.seed(seed)
    nt <- data.frame(TF = paste0(e_tf, "_A"), Gene = paste0(e_gene, "_B"),
                     ctrl = ctrl_w, cond = cond_w, stringsAsFactors = FALSE)
    cp <- nt[nt$ctrl >= 0, 1:3]
    el <- data.frame(red = cp[, 2], blue = cp[, 1], weights = cp[, 3])
    invisible(capture.output({
      co <- netZooR:::createCondorObject(el)
      co <- netZooR:::condorCluster(co, project = FALSE)
    }))
    cm <- c(co$red.memb[, 2], co$blue.memb[, 2])
    names(cm) <- c(as.character(co$red.memb[, 1]), as.character(co$blue.memb[, 1]))
    pt <- nt[nt$ctrl >= 0 & nt$cond >= 0, ]
    pg <- igraph::graph_from_edgelist(as.matrix(pt[, 1:2]), directed = TRUE)
    unc <- setdiff(V(pg)$name, names(cm))
    if (length(unc)) {
      um <- sample(seq_len(max(cm)), length(unc), replace = TRUE); names(um) <- unc; cm <- c(cm, um)
    }
    dwbm <- netZooR:::alpacaComputeDWBMmatmScale(pt, cm[V(pg)$name])
    na <- nrow(dwbm); nb <- ncol(dwbm)
    B <- matrix(0, na + nb, na + nb)
    B[seq_len(na), na + seq_len(nb)] <- dwbm
    B[na + seq_len(nb), seq_len(na)] <- t(dwbm)
    invisible(capture.output(memb <- as.vector(netZooR:::alpacaGenLouvain(B))))
    if (length(memb) != na + nb) stop("alpacaGenLouvain returned no partition")
    names(memb) <- c(rownames(dwbm), colnames(dwbm))
    rm(B)
    mods <- sort(unique(memb))
    mt <- do.call(rbind, lapply(mods, function(i) {
      nm <- names(memb)[memb == i]
      a <- nm[grepl("_A$", nm)]; b <- nm[grepl("_B$", nm)]
      q <- if (length(a) && length(b)) sum(dwbm[a, b, drop = FALSE]) else 0
      data.frame(Module = i, n_tf = length(a), n_gene = length(b), Q = q)
    }))
    mt$n_nodes <- mt$n_tf + mt$n_gene
    mt$density <- ifelse(mt$n_tf > 0 & mt$n_gene > 0, mt$Q / (mt$n_tf * mt$n_gene), 0)
    out <- list(memb = memb, mod = mt, Q_total = sum(mt$Q),
                Q_max = max(mt$Q),
                dens_max = max(c(0, mt$density[mt$n_tf > 0 & mt$n_gene > 0 & mt$n_nodes >= 5])))
    if (keep_dwbm) out$dwbm <- dwbm
    out
  }
  group_means <- function(fxs) {
    ctl <- setdiff(samp, fxs)
    list(ctrl = rowMeans(lioness_mat[, ctl, drop = FALSE]),
         cond = rowMeans(lioness_mat[, fxs, drop = FALSE]))
  }

  # ---- all within-timepoint relabelings ----
  by_tp <- split(samp, tp)
  stopifnot(length(by_tp) == 2, all(lengths(by_tp) == 6),
            all(vapply(by_tp, function(s) sum(cond[match(s, samp)] == "FXS"), 0) == 3))
  cmb <- lapply(by_tp, function(s) combn(s, 3, simplify = FALSE))
  gr  <- expand.grid(a = seq_along(cmb[[1]]), b = seq_along(cmb[[2]]))
  relabel <- lapply(seq_len(nrow(gr)), function(i) sort(c(cmb[[1]][[gr$a[i]]], cmb[[2]][[gr$b[i]]])))
  true_fxs <- sort(samp[cond == "FXS"])
  is_true <- vapply(relabel, identical, logical(1), true_fxs)
  stopifnot(sum(is_true) == 1)
  set.seed(777); order_perm <- sample(which(!is_true))
  N_PERM <- min(env_num("FXS_P5_ALPACA_PERM", 100), length(order_perm))
  cat("Distinct relabelings:", length(relabel), "| running", N_PERM,
      "| smallest attainable p:", signif(1 / (N_PERM + 1), 3), "\n")

  # ---- observed ----
  obs_f <- ckpt("p5_alpaca_obs.rds")
  if (file.exists(obs_f)) { A_obs <- readRDS(obs_f); cat("Observed ALPACA loaded.\n") } else {
    cat("[", ts(), "] Observed ALPACA (control = reference)...\n")
    gm <- group_means(true_fxs)
    t1 <- Sys.time()
    A_obs <- run_alpaca(gm$ctrl, gm$cond, seed = 4000, keep_dwbm = TRUE)
    cat("One ALPACA run took", round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1), "min\n")
    saveRDS(A_obs, obs_f); rm(gm)
  }
  cat("Observed: modules", nrow(A_obs$mod), "| with >=5 nodes", sum(A_obs$mod$n_nodes >= 5),
      "| Q_total", signif(A_obs$Q_total, 4), "\n")

  # ---- permutations ----
  sf <- ckpt("p5_alpaca_state.rds")
  st <- if (file.exists(sf)) readRDS(sf) else list(done = 0, stats = data.frame())
  if (st$done > 0) cat("Resuming after relabeling", st$done, "\n")
  t0 <- Sys.time(); start_done <- st$done
  while (st$done < N_PERM) {
    j <- st$done + 1
    fxs_j <- relabel[[order_perm[j]]]
    gm <- group_means(fxs_j)
    r <- tryCatch(run_alpaca(gm$ctrl, gm$cond, seed = 4000 + j),
                  error = function(e) { cat("  relabeling", j, "failed:", conditionMessage(e), "\n"); NULL })
    if (!is.null(r))
      st$stats <- rbind(st$stats, data.frame(
        perm = j, fxs_samples = paste(fxs_j, collapse = ";"),
        overlap_with_true = length(intersect(fxs_j, true_fxs)),
        Q_total = r$Q_total, Q_max = r$Q_max, dens_max = r$dens_max,
        n_modules_5plus = sum(r$mod$n_nodes >= 5)))
    st$done <- j
    saveRDS(st, sf)
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("[%s] relabeling %d/%d | %.1f min | ETA %.0f min\n", ts(), j, N_PERM,
                el, el / (j - start_done) * (N_PERM - j)))
    rm(gm, r); gc(verbose = FALSE)
  }

  ps <- st$stats; nP <- nrow(ps)
  cat("Successful relabelings:", nP, "\n")
  write.csv(ps, file.path(OUT, "ALPACA_labelperm_null_runs.csv"), row.names = FALSE)

  # ---- module-level, family-wise (max-statistic) p-values ----
  mt <- A_obs$mod
  mt$p_global_maxQ <- vapply(mt$Q, function(q) (1 + sum(ps$Q_max >= q)) / (nP + 1), numeric(1))
  elig <- mt$n_tf > 0 & mt$n_gene > 0 & mt$n_nodes >= 5
  mt$p_global_maxDensity <- NA_real_
  mt$p_global_maxDensity[elig] <- vapply(mt$density[elig],
    function(dn) (1 + sum(ps$dens_max >= dn)) / (nP + 1), numeric(1))
  mt <- mt[order(mt$p_global_maxQ, -mt$Q), ]
  p_total <- (1 + sum(ps$Q_total >= A_obs$Q_total)) / (nP + 1)
  cat("Global differential modularity: observed Q_total =", signif(A_obs$Q_total, 4),
      "| null median", signif(median(ps$Q_total), 4), "| p =", signif(p_total, 3), "\n")
  cat("Top modules (family-wise p):\n"); print(head(mt, 10), row.names = FALSE)
  write.csv(mt, file.path(OUT, "ALPACA_labelperm_modules.csv"), row.names = FALSE)

  # ---- corrected node table (ALPACA standard score + raw contribution) ----
  memb <- A_obs$memb; dw <- A_obs$dwbm
  node_tab <- do.call(rbind, lapply(sort(unique(memb)), function(i) {
    nm <- names(memb)[memb == i]; a <- nm[grepl("_A$", nm)]; b <- nm[grepl("_B$", nm)]
    if (!(length(a) && length(b))) return(data.frame(Node = nm, Module = i, Contribution = 0, Score_norm = 0))
    sub <- dw[a, b, drop = FALSE]; den <- sum(sub)
    contrib <- c(rowSums(sub), colSums(sub))
    data.frame(Node = c(a, b), Module = i, Contribution = contrib,
               Score_norm = if (den != 0) contrib / den else 0)
  }))
  node_tab$Gene <- sub("_[AB]$", "", node_tab$Node)
  node_tab$Type <- ifelse(grepl("_A$", node_tab$Node), "TF", "Target")
  node_tab <- node_tab[order(-node_tab$Contribution), ]
  write.csv(node_tab, file.path(OUT, "ALPACA_corrected_nodes.csv"), row.names = FALSE)

  # Cytoscape tables for the corrected modules (read by 14_cytoscape_local.R):
  # the five modules with the largest Q among those with >=5 nodes, plus any
  # family-wise significant module.
  top_m <- unique(c(head(mt$Module[mt$n_nodes >= 5 & mt$n_tf > 0 & mt$n_gene > 0][
                      order(-mt$Q[mt$n_nodes >= 5 & mt$n_tf > 0 & mt$n_gene > 0])], 5),
                    mt$Module[mt$p_global_maxQ < 0.05]))
  cy_nodes <- node_tab %>% filter(Module %in% top_m) %>%
    group_by(Gene) %>% slice_max(Contribution, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(id = Gene, Module = as.character(Module), Score = round(Score_norm, 6),
              Contribution = Contribution, Type, label = Gene,
              NodeColor = ifelse(Type == "TF", "#E74C3C", "#3498DB"))
  cy_edges <- edge_list_panda %>% filter(TF %in% cy_nodes$id, Target %in% cy_nodes$id) %>%
    transmute(TF, Gene = Target, Score)
  write.csv(cy_nodes, file.path(OUT, "ALPACA_corrected_cyto_nodes.csv"), row.names = FALSE)
  fwrite(cy_edges, file.path(OUT, "ALPACA_corrected_cyto_edges.csv"))
  cat("Cytoscape tables: modules", paste(top_m, collapse = ", "), "|", nrow(cy_nodes), "nodes,",
      nrow(cy_edges), "edges\n")

  prev4 <- c("BARHL1_A", "OTX1_A", "NEFL_B", "MPZL2_B", "TAF4B_B")
  placement <- data.frame(Node = prev4, Module = unname(memb[prev4]))
  placement <- merge(placement, mt[, c("Module", "n_nodes", "Q", "p_global_maxQ")], by = "Module", all.x = TRUE)
  cat("Where the former module-4 genes fall now:\n"); print(placement, row.names = FALSE)
  write.csv(placement, file.path(OUT, "ALPACA_former_module4_placement.csv"), row.names = FALSE)

  # ---- module-local centrality of known TFs, percentile among TFs only ----
  lc <- do.call(rbind, lapply(sort(unique(memb)), function(i) {
    g <- unique(sub("_[AB]$", "", names(memb)[memb == i]))
    tfs_i <- unique(sub("_A$", "", names(memb)[memb == i & grepl("_A$", names(memb))]))
    if (length(g) < 5 || length(tfs_i) < 5) return(NULL)
    se <- edge_list_panda[edge_list_panda$TF %in% tfs_i & edge_list_panda$Target %in% g, ]
    if (nrow(se) < 3) return(NULL)
    s <- tapply(se$Score, se$TF, sum)
    s <- s[names(s) %in% tfs_i]
    full <- setNames(numeric(length(tfs_i)), tfs_i); full[names(s)] <- s
    data.frame(Gene = names(full), Module = i, n_TFs_in_module = length(full),
               Local_Strength = as.numeric(full),
               Pct_among_TFs = rank(as.numeric(full)) / length(full))
  }))
  if (!is.null(lc)) {
    write.csv(lc, file.path(OUT, "ALPACA_local_centrality_TFs.csv"), row.names = FALSE)
    cat("Known TFs — within-module strength percentile among TFs:\n")
    print(lc[lc$Gene %in% known_tf, ], row.names = FALSE)
  }

  # ---- figures ----
  png(file.path(FIGO, "ALPACA_labelperm_null.png"), width = 2000, height = 900, res = 200)
  par(mfrow = c(1, 2))
  hist(ps$Q_total, breaks = 30, col = "#95A5A6", border = "white",
       xlim = range(c(ps$Q_total, A_obs$Q_total)),
       main = sprintf("Total differential modularity (p = %.3f)", p_total),
       xlab = "Q_total under relabeling")
  abline(v = A_obs$Q_total, col = "#E74C3C", lwd = 2)
  hist(ps$Q_max, breaks = 30, col = "#95A5A6", border = "white",
       xlim = range(c(ps$Q_max, mt$Q)),
       main = "Best module under relabeling vs observed modules", xlab = "Module Q")
  abline(v = head(sort(mt$Q, decreasing = TRUE), 5), col = "#E74C3C", lwd = 1.5)
  dev.off()

  alpaca_rows <- data.frame(
    Metric = c("ALPACA relabelings (of 399)", "Observed modules", "Observed modules with >=5 nodes",
               "Global differential modularity p", "Best module family-wise p (Q)",
               "Modules family-wise p<0.05 (Q)", "Best module family-wise p (density, >=5 nodes)"),
    Value = c(nP, nrow(mt), sum(mt$n_nodes >= 5), signif(p_total, 3),
              signif(min(mt$p_global_maxQ), 3), sum(mt$p_global_maxQ < 0.05),
              signif(suppressWarnings(min(mt$p_global_maxDensity, na.rm = TRUE)), 3)))
  write.csv(alpaca_rows, file.path(OUT, "ALPACA_labelperm_summary.csv"), row.names = FALSE)
  print(alpaca_rows, row.names = FALSE)
}

# ============================================================
# STAGE 4 — SUMMARY
# ============================================================
if ("summary" %in% STAGES) {
  cat("\n=== PART 5 SUMMARY ===\n")
  parts <- list(GENIE3 = "GENIE3_matched_null_summary.csv",
                PANDA  = "PANDA_robust_null_summary.csv",
                ALPACA = "ALPACA_labelperm_summary.csv")
  tab <- do.call(rbind, lapply(names(parts), function(n) {
    f <- file.path(OUT, parts[[n]])
    if (file.exists(f)) cbind(Test = n, read.csv(f, stringsAsFactors = FALSE)) else NULL
  }))
  if (!is.null(tab)) {
    print(tab, row.names = FALSE)
    write.csv(tab, file.path(OUT, "part5_summary.csv"), row.names = FALSE)
  }

  # Retire the invalid v5 ALPACA significance outputs
  old_dir <- file.path("results", "ALPACA", "superseded_v5")
  old_f <- c(file.path("results", "ALPACA", c("ALPACA_module_significance.csv",
                                              "ALPACA_results_with_significance.csv",
                                              "ALPACA_module_local_centrality.csv")),
             file.path("figures", "ALPACA", "ALPACA_module_significance.png"))
  old_f <- old_f[file.exists(old_f)]
  if (length(old_f)) {
    dir.create(old_dir, showWarnings = FALSE, recursive = TRUE)
    file.rename(old_f, file.path(old_dir, basename(old_f)))
    cat("Moved superseded v5 ALPACA outputs to", old_dir, "\n")
  }

  # Update the headline table in place
  hf <- file.path("results", "v5_headline_numbers.csv")
  if (file.exists(hf)) {
    h <- read.csv(hf, stringsAsFactors = FALSE); h$Value <- as.character(h$Value)
    new_rows <- robust_headline_rows()
    h <- h[!h$Item %in% c("ALPACA nodes pooled FDR<0.05", "Best ALPACA module FDR",
                          "PANDA edges moved by real expression FDR<0.05", new_rows$Item), ]
    h <- rbind(h, new_rows)
    write.csv(h, hf, row.names = FALSE)
    cat("\nHeadline table updated:\n"); print(h, row.names = FALSE)
  }
  cat("\n=== PART 5 COMPLETE ===\n")
}
