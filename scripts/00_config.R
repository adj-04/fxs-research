# ============================================================
# FXS PIPELINE v5 — shared config, sourced by every part.
#
# Works unchanged on Colab and locally:
#   FXS_ROOT  = where data/, results/, figures/, checkpoints/ live
#               (Colab: a Google Drive folder, so work survives
#               disconnects). Default: current working directory.
#   FXS_REPO  = the cloned repo (for tf_list.txt and scripts).
#               Default: parent of this script's folder, else FXS_ROOT.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || identical(a, "")) b else a

FXS_ROOT <- Sys.getenv("FXS_ROOT", unset = getwd())
FXS_REPO <- Sys.getenv("FXS_REPO", unset = FXS_ROOT)
dir.create(FXS_ROOT, recursive = TRUE, showWarnings = FALSE)
setwd(FXS_ROOT)
cat("FXS_ROOT:", FXS_ROOT, "\nFXS_REPO:", FXS_REPO, "\n")

DATA_DIR        <- "data"
CHECKPOINTS_DIR <- "checkpoints"
NETWORKS_DIR    <- "networks"
LOG_DIR         <- "logs"

RESULTS_ALPACA  <- file.path("results", "ALPACA")
RESULTS_DE      <- file.path("results", "DE")
RESULTS_GENIE3  <- file.path("results", "GENIE3")
RESULTS_LIONESS <- file.path("results", "LIONESS")
RESULTS_METH    <- file.path("results", "methylation")
RESULTS_OTHER   <- file.path("results", "other")
RESULTS_PANDA   <- file.path("results", "PANDA")
RESULTS_VALID   <- file.path("results", "validation")
RESULTS_TISSUE  <- file.path("results", "tissue_generality")

FIG_ALPACA  <- file.path("figures", "ALPACA")
FIG_DE      <- file.path("figures", "DE")
FIG_GENIE3  <- file.path("figures", "GENIE3")
FIG_LIONESS <- file.path("figures", "LIONESS")
FIG_METH    <- file.path("figures", "methylation")
FIG_OTHER   <- file.path("figures", "other")
FIG_PANDA   <- file.path("figures", "PANDA")
FIG_VALID   <- file.path("figures", "validation")
FIG_TISSUE  <- file.path("figures", "tissue_generality")

for (d in c(DATA_DIR, CHECKPOINTS_DIR, NETWORKS_DIR, LOG_DIR,
            RESULTS_ALPACA, RESULTS_DE, RESULTS_GENIE3, RESULTS_LIONESS,
            RESULTS_METH, RESULTS_OTHER, RESULTS_PANDA, RESULTS_VALID,
            RESULTS_TISSUE, FIG_ALPACA, FIG_DE, FIG_GENIE3, FIG_LIONESS,
            FIG_METH, FIG_OTHER, FIG_PANDA, FIG_VALID, FIG_TISSUE))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# ANALYSIS CONFIG — every v5 change is switchable so the v4
# numbers can be reproduced for the paper's sensitivity table.
# ------------------------------------------------------------

env_num  <- function(k, d) as.numeric(Sys.getenv(k, unset = d))
env_bool <- function(k, d) toupper(Sys.getenv(k, unset = d)) %in% c("TRUE", "1", "YES")

# Colab free = 2 vCPU / ~12.7 GB. Forking copies memory, so cap cores.
N_CORES <- min(parallel::detectCores(), env_num("FXS_CORES", 2))

# --- Target universe ---
# v4: GENIE3 used TFs + 13 DEGs (a TF-TF network); PANDA targets were the
#     top 3000 genes by DE p-value (selection on the outcome).
# v5: both use a label-independent universe = top-variance genes (after
#     timepoint adjustment) plus DEGs, so GENIE3 becomes a TF->target map.
TARGET_SELECTION  <- Sys.getenv("FXS_TARGET_SELECTION", unset = "variance") # or "pvalue" (v4)
N_GENIE3_TARGETS  <- env_num("FXS_N_GENIE3_TARGETS", 5000)  # non-TF targets
N_TARGET_POOL     <- env_num("FXS_N_TARGET_POOL",    3000)  # PANDA motif target pool
# Edge-weight quantile for the GENIE3 igraph network (Step 9). v4 kept the
# top decile of ~1.7M edges; with ~8M edges the top 1% keeps a similar,
# tractable size for betweenness and layout.
GENIE3_EDGE_QUANTILE <- env_num("FXS_GENIE3_EDGE_QUANTILE",
                                if (TARGET_SELECTION == "pvalue") 0.90 else 0.99)

# --- Covariates ---
ADJUST_TIMEPOINT_FOR_NETWORKS <- env_bool("FXS_ADJUST_TIMEPOINT", "TRUE")
METH_ADJUST_CELLTYPE          <- env_bool("FXS_METH_ADJUST_CELLTYPE", "TRUE")

# --- Permutation nulls ---
N_PERM_GENIE3  <- env_num("FXS_N_PERM_GENIE3", 20)   # pooled null -> 20 is ample
N_TREES_PERM   <- env_num("FXS_N_TREES_PERM", 250)
N_TREES_OBS    <- env_num("FXS_N_TREES_OBS", 1000)
N_PERM_PANDA   <- env_num("FXS_N_PERM_PANDA", 20)
N_PERM_ALPACA  <- env_num("FXS_N_PERM_ALPACA", 50)
N_BOOT_LIONESS <- env_num("FXS_N_BOOT_LIONESS", 1000)

set.seed(42)
known_tf <- c("SP1", "NRF1", "TFAP2A", "REST", "SIX3", "FOXG1")

# ------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------

ckpt <- function(name) file.path(CHECKPOINTS_DIR, name)

# Pooled null accumulated as a histogram on a fine log grid: memory is
# constant however many null values are added, and p-value resolution is
# 1 / (total null values) instead of 1 / (n_perm + 1). Values <= floor are
# counted in the lowest bin.
#   log = TRUE : grid on log10(x), for positive scores (GENIE3 weights)
#   log = FALSE: linear grid on [lo, hi], for signed scores (PANDA, ALPACA)
make_null_hist <- function(lo, hi, log = FALSE, n_bins = 200000) {
  list(breaks = seq(lo, hi, length.out = n_bins + 1),
       counts = numeric(n_bins), n = 0, lo = lo, hi = hi, log = log)
}
.hist_bin <- function(h, x) {
  v <- if (h$log) log10(pmax(x, 10^h$lo)) else x
  v <- pmin(pmax(v, h$lo), h$hi - 1e-12)
  findInterval(v, h$breaks, rightmost.closed = TRUE, all.inside = TRUE)
}
add_to_null_hist <- function(h, x) {
  x <- x[is.finite(x)]
  h$counts <- h$counts + tabulate(.hist_bin(h, x), nbins = length(h$counts))
  h$n <- h$n + length(x)
  h
}
# Empirical P(null >= x). Conservative: the whole bin containing x counts.
null_hist_p <- function(h, x) {
  tail_counts <- rev(cumsum(rev(h$counts)))
  (tail_counts[.hist_bin(h, x)] + 1) / (h$n + 1)
}
null_hist_quantile <- function(h, q) {
  i <- which(cumsum(h$counts) / h$n >= q)[1]
  b <- h$breaks[i + 1]
  if (h$log) 10^b else b
}

# v5.1: headline rows from the robust nulls (scripts/part5_robust_nulls.R).
# Used by part 4 (fresh runs) and part 5 (updates the existing table).
robust_headline_rows <- function(dir = file.path("results", "robust_nulls")) {
  rd <- function(f) { f <- file.path(dir, f); if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL }
  gv <- function(df, m) if (is.null(df)) NULL else as.character(df$Value[match(m, df$Metric)])
  g <- rd("GENIE3_matched_null_summary.csv"); p <- rd("PANDA_robust_null_summary.csv")
  a <- rd("ALPACA_labelperm_summary.csv")
  rows <- list(
    c("GENIE3 edges FDR<0.05 (matched-tree null)", gv(g, "Observed edges FDR<0.05")),
    c("GENIE3 null calibration: held-out fraction p<0.05", gv(g, "Calibration: null fraction p<0.05 (mean)")),
    c("PANDA observed fraction above prior-null q95 (all edges)", gv(p, "Observed fraction above null q95")),
    c("PANDA r(observed, shuffled expression)", gv(p, "Correlation observed vs shuffled-expression network (mean)")),
    c("PANDA r(shuffled, shuffled)", gv(p, "Correlation between shuffled-expression networks (mean)")),
    c("PANDA edges moved by real expression FDR<0.05", gv(p, "Edges moved by real expression, FDR<0.05")),
    c("ALPACA relabelings", gv(a, "ALPACA relabelings (of 399)")),
    c("ALPACA global differential modularity p", gv(a, "Global differential modularity p")),
    c("Best ALPACA module family-wise p", gv(a, "Best module family-wise p (Q)")))
  rows <- Filter(function(r) length(r) == 2 && !is.na(r[2]), rows)
  if (!length(rows)) return(NULL)
  data.frame(Item = vapply(rows, `[`, "", 1), Value = vapply(rows, `[`, "", 2))
}

log_session <- function(part) {
  writeLines(capture.output(sessionInfo()),
             file.path(LOG_DIR, paste0("sessionInfo_", part, ".txt")))
}

cat("Config loaded | cores:", N_CORES,
    "| target selection:", TARGET_SELECTION,
    "| GENIE3 targets:", N_GENIE3_TARGETS,
    "| timepoint-adjusted networks:", ADJUST_TIMEPOINT_FOR_NETWORKS,
    "| cell-type-adjusted methylation:", METH_ADJUST_CELLTYPE, "\n")
