# ============================================================
# FXS PIPELINE v5.1 — PART 6 : REVIEW CHECKS
#
# Answers a second round of review points. Needs checkpoint_step16.RData
# and, for some stages, Part 5 state files (p5_genie3_state.rds,
# p5_alpaca_state.rds, p5_alpaca_obs.rds). Internet is needed for
# CollecTRI (OmniPath) and ChEA3.
#
# Stages (FXS_PART6_STAGES, comma-separated; default all):
#   tfact        TF activity on donor-blocked limma-voom statistics in the
#                validation cohort, next to the DESeq2 Wald statistics used
#                in Part 4.
#   rest         Negative controls for ENCODE/ChEA3 and CollecTRI
#                self-recovery: (a) top targets from GENIE3 networks built on
#                gene-shuffled expression, (b) random gene sets from the same
#                universe, (c) PANDA's motif prior alone (circularity check).
#   lioness_mde  Minimum detectable effect of the per-edge LIONESS test at
#                n = 6 vs 6 and at larger hypothetical n.
#   genie3_power Spike-in: targets made to co-express with one TF at a set
#                correlation; fraction of spiked edges reaching FDR < 0.05
#                against the Part 5 pooled null.
#   alpaca399    Completes all 399 ALPACA relabelings in the Part 5 state file
#                (then rerun Part 5 with FXS_PART5_STAGES=alpaca,summary and
#                FXS_P5_ALPACA_PERM=399 to regenerate its tables).
#   alpaca_power Spike-in: FXS-only increase of a TF x gene block of LIONESS
#                edges; detection by the same label-permutation test.
#   summary
#
# Parallel: FXS_P6_WORKERS (default 4) forked workers for ALPACA stages.
# All loops checkpoint and resume.
# ============================================================

suppressPackageStartupMessages({
  library(GENIE3); library(netZooR); library(igraph); library(limma)
  library(dplyr); library(ggplot2); library(data.table); library(parallel)
  library(httr); library(jsonlite)
})

source(file.path(Sys.getenv("FXS_REPO", "."), "scripts", "00_config.R"))
load(ckpt("checkpoint_step16.RData"))
source(file.path(Sys.getenv("FXS_REPO", "."), "scripts", "00_config.R"))
log_session("part6")
cat("Checkpoint loaded.\n")

STAGES <- trimws(strsplit(Sys.getenv("FXS_PART6_STAGES",
  "tfact,rest,lioness_mde,genie3_power,alpaca399,alpaca_power,summary"), ",")[[1]])
W <- max(1, env_num("FXS_P6_WORKERS", 4))
OUT  <- file.path("results", "review_checks")
FIGO <- file.path("figures", "review_checks")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGO, recursive = TRUE, showWarnings = FALSE)
cat("Stages:", paste(STAGES, collapse = ", "), "| workers:", W, "\n")
ts <- function() format(Sys.time(), "%H:%M:%S")

get_collectri <- function() {
  f <- ckpt("p6_collectri.rds")
  if (file.exists(f)) return(readRDS(f))
  suppressPackageStartupMessages(library(decoupleR))
  ct <- as.data.frame(decoupleR::get_collectri(organism = "human", split_complexes = FALSE))
  saveRDS(ct, f); ct
}

# ============================================================
# STAGE tfact — TF activity, donor-blocked
# ============================================================
if ("tfact" %in% STAGES) {
  cat("\n=== TF activity on donor-blocked statistics ===\n")
  suppressPackageStartupMessages(library(decoupleR))
  collectri <- get_collectri()
  ulm <- function(stat, sym) {
    d <- data.frame(Symbol = sym, stat = stat) %>%
      filter(!is.na(stat), !is.na(Symbol), Symbol != "") %>%
      group_by(Symbol) %>% slice_max(abs(stat), n = 1, with_ties = FALSE) %>% ungroup()
    m <- matrix(d$stat, ncol = 1, dimnames = list(d$Symbol, "stat"))
    r <- as.data.frame(decoupleR::run_ulm(mat = m, network = collectri, .source = "source",
                                          .target = "target", .mor = "mor", minsize = 5))
    r[, c("source", "score", "p_value")]
  }
  a_org <- ulm(res_df$stat, res_df$Symbol)
  a_wald <- ulm(deg_valid$stat, deg_valid$Symbol)
  a_blk <- ulm(v_tt$t, v_tt$Symbol)
  nm <- function(a, s) setNames(a, c("source", paste0("score_", s), paste0("p_", s)))
  ta <- Reduce(function(x, y) merge(x, y, by = "source", all = TRUE),
               list(nm(a_org, "organoid"), nm(a_wald, "valid_DESeq2"),
                    nm(a_blk, "valid_voom_donorblock")))
  for (s in c("organoid", "valid_DESeq2", "valid_voom_donorblock"))
    ta[[paste0("FDR_", s)]] <- p.adjust(ta[[paste0("p_", s)]], "BH")
  ta$Known_FXS <- ta$source %in% known_tf
  ta <- ta[order(ta$p_valid_voom_donorblock), ]
  write.csv(ta, file.path(OUT, "TF_activity_donorblock.csv"), row.names = FALSE)

  focus <- c(known_tf, "CTNNB1")
  cat("Reference TFs and CTNNB1:\n")
  print(ta[ta$source %in% focus, c("source", "score_organoid", "FDR_organoid",
          "score_valid_DESeq2", "FDR_valid_DESeq2",
          "score_valid_voom_donorblock", "FDR_valid_voom_donorblock")], row.names = FALSE)

  both <- function(col) ta$source[which(ta$FDR_organoid < 0.05 & ta[[col]] < 0.05 &
                                          sign(ta$score_organoid) == sign(ta[[sub("FDR_", "score_", col)]]))]
  rho <- function(s) suppressWarnings(cor(ta$score_organoid, ta[[paste0("score_", s)]],
                                          method = "spearman", use = "complete.obs"))
  tf_sum <- data.frame(
    Metric = c("Donor correlation (duplicateCorrelation)",
               "TFs FDR<0.05, validation DESeq2 (Wald)",
               "TFs FDR<0.05, validation limma-voom (donor-blocked)",
               "Same-direction FDR<0.05 in both cohorts, DESeq2",
               "Same-direction FDR<0.05 in both cohorts, donor-blocked",
               "Spearman organoid vs validation, DESeq2",
               "Spearman organoid vs validation, donor-blocked"),
    Value = c(round(v_corr, 4),
              sum(ta$FDR_valid_DESeq2 < 0.05, na.rm = TRUE),
              sum(ta$FDR_valid_voom_donorblock < 0.05, na.rm = TRUE),
              paste(both("FDR_valid_DESeq2"), collapse = ";"),
              paste(both("FDR_valid_voom_donorblock"), collapse = ";"),
              round(rho("valid_DESeq2"), 4), round(rho("valid_voom_donorblock"), 4)))
  print(tf_sum, row.names = FALSE)
  write.csv(tf_sum, file.path(OUT, "TF_activity_donorblock_summary.csv"), row.names = FALSE)
}

# ============================================================
# STAGE rest — negative controls for regulon self-recovery
# ============================================================
if ("rest" %in% STAGES) {
  cat("\n=== Self-recovery negative controls (ChEA3 ENCODE, CollecTRI) ===\n")
  chea3_query <- function(gene_set, query_name = "query") {
    gene_set <- unique(toupper(as.character(gene_set)))
    gene_set <- gene_set[!is.na(gene_set) & gene_set != ""]
    if (length(gene_set) < 5) return(NULL)
    for (attempt in 1:3) {
      resp <- tryCatch(POST("https://maayanlab.cloud/chea3/api/enrich/",
                            body = toJSON(list(query_name = query_name, gene_set = gene_set),
                                          auto_unbox = TRUE),
                            encode = "raw", content_type_json(), timeout(120)),
                       error = function(e) NULL)
      if (!is.null(resp) && status_code(resp) == 200) break
      Sys.sleep(5)
    }
    if (is.null(resp) || status_code(resp) != 200) return(NULL)
    res <- tryCatch(fromJSON(content(resp, as = "text", encoding = "UTF-8")), error = function(e) NULL)
    enc <- res[["ENCODE--ChIP-seq"]]
    if (is.null(enc) || !"TF" %in% colnames(enc)) return(NULL)
    toupper(as.character(enc$TF))                      # ranked ENCODE TFs
  }
  enc_rank <- function(ranked, tf) { r <- which(ranked == toupper(tf)); if (length(r)) r[1] else NA_integer_ }

  obs_sr <- read.csv(file.path(RESULTS_VALID, "ChEA3_ENCODE_selfrecovery.csv"), stringsAsFactors = FALSE)
  tfs_g <- unique(obs_sr$TF[obs_sr$Source == "GENIE3"])
  tfs_p <- unique(obs_sr$TF[obs_sr$Source == "PANDA"])
  cat("TFs:", paste(unique(obs_sr$TF), collapse = ", "), "\n")

  sf <- ckpt("p6_rest_state.rds")
  st <- if (file.exists(sf)) readRDS(sf) else list(perm = list(), rand_g = list(), rand_p = list(), prior = NULL)

  # (a) GENIE3 on gene-shuffled expression (same seeds as Part 5)
  N_G <- env_num("FXS_P6_GENIE3_PERM", 5)
  permute_rows <- function(m) { o <- t(apply(m, 1, sample)); dimnames(o) <- dimnames(m); o }
  for (p in seq_len(N_G)) {
    if (length(st$perm) >= p && !is.null(st$perm[[p]])) next
    cat(sprintf("[%s] shuffled GENIE3 %d/%d\n", ts(), p, N_G))
    set.seed(5000 + p)
    w <- GENIE3(permute_rows(expr_subset), regulators = tf_in_data,
                nTrees = N_TREES_OBS, nCores = N_CORES, verbose = FALSE)
    tg <- lapply(setNames(tfs_g, tfs_g), function(tf) {
      v <- w[tf, ]; v <- v[names(v) != tf]; names(sort(v, decreasing = TRUE))[1:500] })
    rk <- vapply(tfs_g, function(tf) enc_rank(chea3_query(tg[[tf]][1:200], paste0(tf, "_perm", p)), tf), 1L)
    st$perm[[p]] <- list(targets = tg, enc_rank = rk)
    saveRDS(st, sf); rm(w); gc(verbose = FALSE)
  }

  # (b) random 200-gene sets: one query ranks every ENCODE TF at once
  N_R <- env_num("FXS_P6_RANDOM_SETS", 30)
  uni_g <- colnames(genie3_weights)
  uni_p <- unique(edge_list_panda$Target)
  for (k in seq_len(N_R)) {
    if (length(st$rand_g) < k) {
      set.seed(9000 + k); st$rand_g[k] <- list(chea3_query(sample(uni_g, 200), paste0("randG", k)))
    }
    if (length(st$rand_p) < k) {
      set.seed(9500 + k); st$rand_p[k] <- list(chea3_query(sample(uni_p, 200), paste0("randP", k)))
    }
    saveRDS(st, sf)
    if (k %% 10 == 0) cat(sprintf("[%s] random sets %d/%d\n", ts(), k, N_R))
  }

  # (c) PANDA's motif prior alone
  if (is.null(st$prior)) {
    st$prior <- lapply(setNames(tfs_p, tfs_p), function(tf) {
      tg <- motif_edgelist %>% filter(TF == tf) %>% arrange(desc(Score)) %>% pull(Gene) %>% unique()
      list(targets = head(tg, 500),
           enc_rank = if (length(tg) >= 5) enc_rank(chea3_query(head(tg, 200), paste0(tf, "_prior")), tf) else NA_integer_)
    })
    saveRDS(st, sf)
  }

  rand_rank <- function(lst, tf) vapply(lst, function(r) if (is.null(r)) NA_integer_ else enc_rank(r, tf), 1L)
  emp_p <- function(obs, null) { null <- null[!is.na(null)]
    if (is.na(obs) || !length(null)) NA_real_ else (1 + sum(null <= obs)) / (length(null) + 1) }

  rows <- list()
  for (i in seq_len(nrow(obs_sr))) {
    tf <- obs_sr$TF[i]; src <- obs_sr$Source[i]; o <- obs_sr$Self_rank_ENCODE[i]
    rr <- rand_rank(if (src == "GENIE3") st$rand_g else st$rand_p, tf)
    pr <- if (src == "GENIE3") vapply(st$perm, function(x) x$enc_rank[[tf]], 1L) else NA_integer_
    pri <- if (src == "PANDA" && !is.null(st$prior[[tf]])) st$prior[[tf]]$enc_rank else NA_integer_
    rows[[i]] <- data.frame(TF = tf, Source = src, Observed_ENCODE_rank = o,
      Random_sets_median_rank = median(rr, na.rm = TRUE),
      Random_sets_frac_rank_le_obs = mean(rr <= o, na.rm = TRUE),
      Random_sets_p = emp_p(o, rr),
      Shuffled_GENIE3_ranks = if (src == "GENIE3") paste(pr, collapse = ";") else NA,
      Motif_prior_only_rank = pri)
  }
  ctrl <- do.call(rbind, rows)
  print(ctrl, row.names = FALSE)
  write.csv(ctrl, file.path(OUT, "ChEA3_selfrecovery_negative_controls.csv"), row.names = FALSE)

  # CollecTRI enrichment of the same target sets (top 500, as in Part 4)
  collectri <- get_collectri()
  universe_genes <- unique(c(edge_list_genie3$Target, edge_list_panda$Target,
                             edge_list_genie3$TF, edge_list_panda$TF))
  enr <- function(tf, pred) {
    cur <- intersect(unique(collectri$target[collectri$source == tf]), universe_genes)
    pred <- intersect(unique(pred), universe_genes)
    if (length(cur) < 3 || length(pred) < 3) return(c(fold = NA, p = NA))
    ov <- length(intersect(pred, cur)); ex <- length(pred) * length(cur) / length(universe_genes)
    c(fold = ov / ex, p = phyper(ov - 1, length(cur), length(universe_genes) - length(cur),
                                 length(pred), lower.tail = FALSE))
  }
  get_obs <- function(tf, src) {
    if (src == "GENIE3") edge_list_genie3 %>% filter(TF == tf) %>% arrange(desc(Weight)) %>% head(500) %>% pull(Target)
    else edge_list_panda %>% filter(TF == tf) %>% arrange(desc(Score)) %>% head(500) %>% pull(Target)
  }
  ce <- do.call(rbind, lapply(intersect(unique(obs_sr$TF), unique(collectri$source)), function(tf) {
    do.call(rbind, lapply(c("GENIE3", "PANDA"), function(src) {
      o <- enr(tf, as.character(get_obs(tf, src)))
      sh <- if (src == "GENIE3") sapply(st$perm, function(x) enr(tf, x$targets[[tf]])["p"]) else NA
      pri <- if (src == "PANDA" && !is.null(st$prior[[tf]])) enr(tf, st$prior[[tf]]$targets)["p"] else NA
      data.frame(TF = tf, Source = src, Observed_fold = o[["fold"]], Observed_p = o[["p"]],
                 Shuffled_GENIE3_min_p = if (all(is.na(sh))) NA else min(sh, na.rm = TRUE),
                 Motif_prior_only_p = pri)
    }))
  }))
  print(ce, row.names = FALSE)
  write.csv(ce, file.path(OUT, "CollecTRI_selfrecovery_negative_controls.csv"), row.names = FALSE)
}

# ============================================================
# STAGE lioness_mde — minimum detectable effect, per-edge LIONESS test
# ============================================================
if ("lioness_mde" %in% STAGES) {
  cat("\n=== LIONESS minimum detectable effect ===\n")
  samp <- colnames(lioness_mat)
  des <- model.matrix(~ timepoint[samp] + condition[samp])
  colnames(des) <- c("Intercept", "timepointD105", "conditionFXS")
  fit <- eBayes(lmFit(lioness_mat, des))
  se  <- sqrt(fit$s2.post) * fit$stdev.unscaled[, "conditionFXS"]
  dft <- fit$df.total
  mu  <- rowMeans(lioness_mat)
  sde <- sqrt(rowSums((lioness_mat - mu)^2) / (ncol(lioness_mat) - 1))
  obs_abs <- abs(fit$coefficients[, "conditionFXS"])
  rm(fit); gc(verbose = FALSE)
  m <- nrow(lioness_mat)
  cat("Edges:", m, "| median residual df (moderated):", round(median(dft), 1), "\n")

  # Effect needed for 80% power at the BH threshold reached when k of m edges
  # are true effects of the same size (threshold ~ 0.05 k / m).
  mde_tab <- do.call(rbind, lapply(c(1, 100, 10000), function(k) {
    pstar <- 0.05 * k / m
    tcrit <- qt(1 - pstar / 2, dft) + qt(0.8, dft)
    mde <- tcrit * se
    data.frame(True_edges_k = k, p_threshold = signif(pstar, 3),
               Median_MDE = round(median(mde), 3),
               Median_MDE_in_edge_SDs = round(median(mde / sde), 2),
               Frac_edges_whose_obs_diff_exceeds_MDE = signif(mean(obs_abs >= mde), 3))
  }))
  print(mde_tab, row.names = FALSE)

  # Same calculation for larger balanced designs (n per group), assuming
  # the same per-edge variance: SE scales with sqrt(6/n), df = 2n - 3.
  n_tab <- do.call(rbind, lapply(c(6, 12, 24, 50, 100), function(n) {
    pstar <- 0.05 * 100 / m; d <- 2 * n - 3
    tcrit <- qt(1 - pstar / 2, d) + qt(0.8, d)
    data.frame(n_per_group = n, Median_MDE_in_edge_SDs_k100 = round(median(tcrit * se * sqrt(6 / n) / sde), 2))
  }))
  print(n_tab, row.names = FALSE)
  write.csv(mde_tab, file.path(OUT, "LIONESS_MDE.csv"), row.names = FALSE)
  write.csv(n_tab, file.path(OUT, "LIONESS_MDE_by_n.csv"), row.names = FALSE)
  rm(se, dft, sde, obs_abs, mu); gc(verbose = FALSE)
}

# ============================================================
# STAGE genie3_power — co-expression spike-in
# ============================================================
if ("genie3_power" %in% STAGES) {
  cat("\n=== GENIE3 spike-in power ===\n")
  g3s <- readRDS(ckpt("p5_genie3_state.rds"))
  h0 <- make_null_hist(lo = -12, hi = 0, log = TRUE)
  hG <- h0; hG$counts <- Reduce(`+`, g3s$counts); hG$n <- sum(unlist(g3s$n))
  obs_p <- null_hist_p(hG, edge_list_genie3$Weight)
  obs_id <- paste(edge_list_genie3$TF, edge_list_genie3$Target, sep = "__")
  m <- length(obs_p)

  rs <- as.numeric(strsplit(Sys.getenv("FXS_P6_SPIKE_R", "0.5,0.7,0.8,0.9,0.95,0.99"), ",")[[1]])
  NT <- env_num("FXS_P6_SPIKE_TARGETS", 50)
  non_tf <- setdiff(rownames(expr_subset), tf_in_data)
  sf <- ckpt("p6_genie3_power.rds")
  res <- if (file.exists(sf)) readRDS(sf) else list()
  for (r in rs) {
    key <- as.character(r)
    if (!is.null(res[[key]])) next
    set.seed(round(7000 + 1000 * r))
    tg <- sample(non_tf, NT); tfs <- sample(tf_in_data, NT, replace = TRUE)
    ex <- expr_subset
    for (i in seq_len(NT)) {
      z <- as.numeric(scale(ex[tfs[i], ])); e <- rnorm(ncol(ex))
      e <- as.numeric(scale(resid(lm(e ~ z))))
      v <- r * z + sqrt(1 - r^2) * e
      ex[tg[i], ] <- mean(ex[tg[i], ]) + sd(ex[tg[i], ]) * v
    }
    w <- GENIE3(ex, regulators = tf_in_data, targets = tg, nTrees = N_TREES_OBS,
                nCores = N_CORES, verbose = FALSE)
    wt <- w[cbind(match(tfs, rownames(w)), match(tg, colnames(w)))]
    p_sp <- null_hist_p(hG, wt)
    # BH across all edges, with the spiked edges replacing their originals
    pp <- obs_p; idx <- match(paste(tfs, tg, sep = "__"), obs_id)
    pp[idx[!is.na(idx)]] <- p_sp[!is.na(idx)]
    if (any(is.na(idx))) pp <- c(pp, p_sp[is.na(idx)])
    fdr_all <- p.adjust(pp, "BH")
    fdr_sp <- c(fdr_all[idx[!is.na(idx)]], if (any(is.na(idx))) tail(fdr_all, sum(is.na(idx))))
    rank_tf <- vapply(seq_len(NT), function(i) rank(-w[, tg[i]])[tfs[i]], numeric(1))
    res[[key]] <- data.frame(r = r, target = tg, tf = tfs, weight = wt, p = p_sp,
                             FDR = fdr_sp, rank_of_true_TF = rank_tf)
    saveRDS(res, sf)
    cat(sprintf("[%s] r = %.2f | median weight %.3f | FDR<0.05: %d/%d | true TF ranked 1st: %d/%d\n",
                ts(), r, median(wt), sum(fdr_sp < 0.05), NT, sum(rank_tf == 1), NT))
    rm(ex, w); gc(verbose = FALSE)
  }
  pw <- do.call(rbind, res)
  pw_sum <- pw %>% group_by(r) %>% summarise(n = n(), median_weight = median(weight),
    detected_FDR05 = mean(FDR < 0.05), true_TF_rank1 = mean(rank_of_true_TF == 1),
    true_TF_top5 = mean(rank_of_true_TF <= 5), .groups = "drop")
  print(as.data.frame(pw_sum), row.names = FALSE)
  write.csv(pw, file.path(OUT, "GENIE3_spikein_edges.csv"), row.names = FALSE)
  write.csv(pw_sum, file.path(OUT, "GENIE3_spikein_power.csv"), row.names = FALSE)
  rm(obs_p, obs_id); gc(verbose = FALSE)
}

# ============================================================
# ALPACA helpers (shared by alpaca399 and alpaca_power)
# ============================================================
if (any(c("alpaca399", "alpaca_power") %in% STAGES)) {
  e_tfA   <- paste0(sub("__.*", "", rownames(lioness_mat)), "_A")
  e_geneB <- paste0(sub(".*__", "", rownames(lioness_mat)), "_B")
  samp <- colnames(lioness_mat)
  tp   <- as.character(timepoint[samp]); cond <- as.character(condition[samp])

  run_alpaca <- function(ctrl_w, cond_w, seed) {
    set.seed(seed)
    nt <- data.frame(TF = e_tfA, Gene = e_geneB, ctrl = ctrl_w, cond = cond_w, stringsAsFactors = FALSE)
    cp <- nt[nt$ctrl >= 0, 1:3]
    el <- data.frame(red = cp[, 2], blue = cp[, 1], weights = cp[, 3]); rm(cp)
    invisible(capture.output({
      co <- netZooR:::createCondorObject(el)
      co <- netZooR:::condorCluster(co, project = FALSE)
    }))
    cm <- c(co$red.memb[, 2], co$blue.memb[, 2])
    names(cm) <- c(as.character(co$red.memb[, 1]), as.character(co$blue.memb[, 1]))
    pt <- nt[nt$ctrl >= 0 & nt$cond >= 0, ]; rm(nt, el, co)
    pg <- igraph::graph_from_edgelist(as.matrix(pt[, 1:2]), directed = TRUE)
    unc <- setdiff(V(pg)$name, names(cm))
    if (length(unc)) { um <- sample(seq_len(max(cm)), length(unc), replace = TRUE); names(um) <- unc; cm <- c(cm, um) }
    dwbm <- netZooR:::alpacaComputeDWBMmatmScale(pt, cm[V(pg)$name])
    na <- nrow(dwbm); nb <- ncol(dwbm)
    B <- matrix(0, na + nb, na + nb)
    B[seq_len(na), na + seq_len(nb)] <- dwbm; B[na + seq_len(nb), seq_len(na)] <- t(dwbm)
    invisible(capture.output(memb <- as.vector(netZooR:::alpacaGenLouvain(B))))
    if (length(memb) != na + nb) stop("alpacaGenLouvain returned no partition")
    names(memb) <- c(rownames(dwbm), colnames(dwbm)); rm(B)
    mt <- do.call(rbind, lapply(sort(unique(memb)), function(i) {
      nm <- names(memb)[memb == i]; a <- nm[grepl("_A$", nm)]; b <- nm[grepl("_B$", nm)]
      data.frame(Module = i, n_tf = length(a), n_gene = length(b),
                 Q = if (length(a) && length(b)) sum(dwbm[a, b, drop = FALSE]) else 0)
    }))
    mt$n_nodes <- mt$n_tf + mt$n_gene
    mt$density <- ifelse(mt$n_tf > 0 & mt$n_gene > 0, mt$Q / (mt$n_tf * mt$n_gene), 0)
    list(memb = memb, mod = mt, Q_total = sum(mt$Q), Q_max = max(mt$Q),
         dens_max = max(c(0, mt$density[mt$n_tf > 0 & mt$n_gene > 0 & mt$n_nodes >= 5])))
  }
  group_means <- function(M, fxs) {
    ctl <- setdiff(samp, fxs)
    list(ctrl = rowMeans(M[, ctl, drop = FALSE]), cond = rowMeans(M[, fxs, drop = FALSE]))
  }
  by_tp <- split(samp, tp)
  cmb <- lapply(by_tp, function(s) combn(s, 3, simplify = FALSE))
  gr  <- expand.grid(a = seq_along(cmb[[1]]), b = seq_along(cmb[[2]]))
  relabel <- lapply(seq_len(nrow(gr)), function(i) sort(c(cmb[[1]][[gr$a[i]]], cmb[[2]][[gr$b[i]]])))
  true_fxs <- sort(samp[cond == "FXS"])
  is_true <- vapply(relabel, identical, logical(1), true_fxs)
  stopifnot(sum(is_true) == 1)
  set.seed(777); order_perm <- sample(which(!is_true))   # identical to Part 5

  # Forked workers share the parent's memory only until R's garbage collector
  # touches it, after which every worker holds its own copy. Drop everything
  # the ALPACA stages do not need before forking (the full checkpoint is
  # several GB).
  keep <- c("lioness_mat", "e_tfA", "e_geneB", "samp", "tp", "cond", "relabel", "true_fxs",
            "order_perm", "run_alpaca", "group_means", "STAGES", "W", "OUT", "FIGO", "ts",
            "ckpt", "CHECKPOINTS_DIR", "env_num", "env_bool", "N_CORES")
  big <- vapply(ls(envir = .GlobalEnv), function(x)
    as.numeric(object.size(get(x, envir = .GlobalEnv))), numeric(1))
  drop <- setdiff(names(big)[big > 20e6], keep)
  cat("Freeing before fork:", paste(drop, collapse = ", "), "\n")
  rm(list = drop, envir = .GlobalEnv); gc(verbose = FALSE)
  cat(sprintf("Remaining R objects: %.1f GB\n", sum(big[setdiff(names(big), drop)]) / 1e9))

  # run relabelings j in parallel on matrix M; returns a stats data.frame
  run_relabelings <- function(M, js, seed_base) {
    out <- mclapply(js, function(j) {
      fx <- relabel[[order_perm[j]]]; gm <- group_means(M, fx)
      r <- tryCatch(run_alpaca(gm$ctrl, gm$cond, seed = seed_base + j), error = function(e) NULL)
      if (is.null(r)) return(NULL)
      data.frame(perm = j, fxs_samples = paste(fx, collapse = ";"),
                 overlap_with_true = length(intersect(fx, true_fxs)),
                 Q_total = r$Q_total, Q_max = r$Q_max, dens_max = r$dens_max,
                 n_modules_5plus = sum(r$mod$n_nodes >= 5))
    }, mc.cores = W, mc.preschedule = FALSE)
    do.call(rbind, out)
  }
}

# ============================================================
# STAGE alpaca399 — complete the exact label-permutation test
# ============================================================
if ("alpaca399" %in% STAGES) {
  cat("\n=== ALPACA: completing all 399 relabelings ===\n")
  sf <- ckpt("p5_alpaca_state.rds")
  st <- readRDS(sf)
  cat("Part 5 state: relabelings done =", st$done, "\n")
  t0 <- Sys.time(); start <- st$done
  while (st$done < length(order_perm)) {
    js <- seq(st$done + 1, min(st$done + W * 3, length(order_perm)))
    st$stats <- rbind(st$stats, run_relabelings(lioness_mat, js, seed_base = 4000))
    st$done <- max(js); saveRDS(st, sf)
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("[%s] relabeling %d/%d | %.1f min | ETA %.0f min\n", ts(), st$done,
                length(order_perm), el, el / (st$done - start) * (length(order_perm) - st$done)))
  }
  cat("Done. Now rerun Part 5 with FXS_PART5_STAGES=alpaca,summary FXS_P5_ALPACA_PERM=399\n")
}

# ============================================================
# STAGE alpaca_power — FXS-only spike-in of a TF x gene block
# ============================================================
if ("alpaca_power" %in% STAGES) {
  cat("\n=== ALPACA spike-in power ===\n")
  deltas <- as.numeric(strsplit(Sys.getenv("FXS_P6_ALPACA_DELTAS", "0.25,0.5,1,2"), ",")[[1]])
  REPS <- env_num("FXS_P6_ALPACA_REPS", 3)
  NREL <- env_num("FXS_P6_ALPACA_NREL", 40)
  NTF <- env_num("FXS_P6_BLOCK_TFS", 15); NG <- env_num("FXS_P6_BLOCK_GENES", 60)
  tf_u <- unique(sub("_A$", "", e_tfA)); g_u <- unique(sub("_B$", "", e_geneB))
  sd_ctrl <- sd(rowMeans(lioness_mat[, setdiff(samp, true_fxs)]))
  cat("Block:", NTF, "TFs x", NG, "genes | delta in units of control edge-score SD =",
      round(sd_ctrl, 3), "| relabelings per scenario:", NREL, "\n")
  row_id <- setNames(seq_len(nrow(lioness_mat)), rownames(lioness_mat))

  sf <- ckpt("p6_alpaca_power.rds")
  res <- if (file.exists(sf)) readRDS(sf) else list()
  for (d in deltas) for (rep in seq_len(REPS)) {
    key <- if (NTF == 15 && NG == 60) sprintf("d%.2f_r%d", d, rep) else sprintf("d%.2f_r%d_b%dx%d", d, rep, NTF, NG)
    if (!is.null(res[[key]])) next
    set.seed(round(8000 + 100 * d + rep))
    btf <- sample(tf_u, NTF); bg <- sample(setdiff(g_u, btf), NG)
    rows <- row_id[paste(rep(btf, each = NG), rep(bg, times = NTF), sep = "__")]
    rows <- rows[!is.na(rows)]
    M <- lioness_mat
    M[rows, true_fxs] <- M[rows, true_fxs] + d * sd_ctrl
    gm <- group_means(M, true_fxs)
    ob <- run_alpaca(gm$ctrl, gm$cond, seed = 4000)
    best <- ob$mod$Module[which.max(ob$mod$Q)]
    best_nodes <- names(ob$memb)[ob$memb == best]
    block_nodes <- c(paste0(btf, "_A"), paste0(bg, "_B"))
    rec <- mean(block_nodes %in% best_nodes)
    nul <- run_relabelings(M, seq_len(NREL), seed_base = 4000)
    res[[key]] <- data.frame(delta = d, rep = rep, block_tfs = NTF, block_genes = NG,
      block_edges = length(rows),
      Q_total = ob$Q_total, Q_max = ob$Q_max,
      p_global = (1 + sum(nul$Q_total >= ob$Q_total)) / (nrow(nul) + 1),
      p_maxQ = (1 + sum(nul$Q_max >= ob$Q_max)) / (nrow(nul) + 1),
      block_in_best_module = rec, n_null = nrow(nul))
    saveRDS(res, sf)
    print(res[[key]], row.names = FALSE)
    rm(M, gm, ob, nul); gc(verbose = FALSE)
  }
  ap <- do.call(rbind, lapply(res, function(x) {
    if (is.null(x$block_tfs)) x <- cbind(x[, 1:2], block_tfs = 15, block_genes = 60, x[, -(1:2)])
    x }))
  ap$frac_edges_spiked <- signif(ap$block_edges / nrow(lioness_mat), 3)
  ap_sum <- ap %>% group_by(block_edges, frac_edges_spiked, delta) %>% summarise(reps = n(), detected_maxQ = mean(p_maxQ < 0.05),
    detected_global = mean(p_global < 0.05), median_block_recovery = median(block_in_best_module),
    .groups = "drop")
  print(as.data.frame(ap_sum), row.names = FALSE)
  write.csv(ap, file.path(OUT, "ALPACA_spikein_runs.csv"), row.names = FALSE)
  write.csv(ap_sum, file.path(OUT, "ALPACA_spikein_power.csv"), row.names = FALSE)
}

# ============================================================
# SUMMARY
# ============================================================
if ("summary" %in% STAGES) {
  cat("\n=== PART 6 OUTPUTS ===\n")
  for (f in list.files(OUT, full.names = TRUE)) cat(" ", f, "\n")
  cat("\n=== PART 6 COMPLETE ===\n")
}
