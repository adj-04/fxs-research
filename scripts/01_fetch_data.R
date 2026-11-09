# ============================================================
# FXS PIPELINE v5 — STEP 0: fetch every input from its source
#
# Nothing has to be downloaded by hand and uploaded. Files already
# present are skipped, so on Colab (FXS_ROOT on Drive) this only
# downloads once. A manifest with sizes and md5 sums is written to
# data/MANIFEST.csv for the paper's Data Availability section.
# ============================================================

source(file.path(Sys.getenv("FXS_REPO", "."), "scripts", "00_config.R"))
options(timeout = 3600)

fetch <- function(url, dest, min_bytes = 1000) {
  if (file.exists(dest) && file.size(dest) >= min_bytes) {
    cat("  have ", basename(dest), " (", round(file.size(dest) / 1e6, 1), " MB)\n", sep = "")
    return(invisible(dest))
  }
  cat("  fetching", basename(dest), "...\n")
  tmp <- paste0(dest, ".part")
  for (attempt in 1:3) {
    ok <- tryCatch({
      download.file(url, tmp, mode = "wb", quiet = TRUE, method = "libcurl"); TRUE
    }, error = function(e) { cat("    attempt", attempt, "failed:", conditionMessage(e), "\n"); FALSE })
    if (ok && file.exists(tmp) && file.size(tmp) >= min_bytes) break
    Sys.sleep(5 * attempt)
  }
  if (!file.exists(tmp) || file.size(tmp) < min_bytes)
    stop("Could not download ", url)
  file.rename(tmp, dest)
  cat("    ok (", round(file.size(dest) / 1e6, 1), " MB)\n", sep = "")
  invisible(dest)
}

geo_counts_url <- function(acc, file)
  paste0("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=",
         acc, "&format=file&file=", file)

cat("\n=== NCBI-reprocessed RNA-seq counts (GRCh38.p13) ===\n")
fetch(geo_counts_url("GSE146339", "GSE146339_raw_counts_GRCh38.p13_NCBI.tsv.gz"),
      file.path(DATA_DIR, "GSE146339_raw_counts_GRCh38.p13_NCBI.tsv.gz"))
fetch(geo_counts_url("GSE218676", "GSE218676_raw_counts_GRCh38.p13_NCBI.tsv.gz"),
      file.path(DATA_DIR, "GSE218676_raw_counts_GRCh38.p13_NCBI.tsv.gz"))
fetch(geo_counts_url("GSE146339", "Human.GRCh38.p13.annot.tsv.gz"),
      file.path(DATA_DIR, "Human.GRCh38.p13.annot.tsv.gz"))

cat("\n=== Sample metadata (shipped in repo, GEO fallback) ===\n")
for (f in c("GSE218676_metadata.csv", "GSE41273_metadata.csv", "tf_list.txt")) {
  src <- file.path(FXS_REPO, "data", f); dst <- file.path(DATA_DIR, f)
  if (!file.exists(dst) && file.exists(src)) file.copy(src, dst)
  cat(" ", f, if (file.exists(dst)) "ok" else "MISSING", "\n")
}

cat("\n=== GSE41273 450K beta matrix (GEO series matrix) ===\n")
meth_rds <- file.path(DATA_DIR, "GSE41273_methylation_matrix.rds")
if (!file.exists(meth_rds)) {
  suppressPackageStartupMessages(library(GEOquery))
  gse <- getGEO("GSE41273", GSEMatrix = TRUE, getGPL = FALSE, destdir = DATA_DIR)
  es  <- gse[[1]]
  beta <- Biobase::exprs(es)
  storage.mode(beta) <- "double"
  rng <- range(beta, na.rm = TRUE)
  cat("  dims:", dim(beta), "| value range:", signif(rng, 3), "\n")
  if (rng[1] < -0.01 || rng[2] > 1.01)
    stop("Series matrix values are not beta values (range outside [0,1]).")
  saveRDS(beta, meth_rds, compress = FALSE)
  pd <- Biobase::pData(es)
  if (!file.exists(file.path(DATA_DIR, "GSE41273_metadata.csv")))
    write.csv(pd, file.path(DATA_DIR, "GSE41273_metadata.csv"))
  rm(gse, es, beta); gc(verbose = FALSE)
} else cat("  have GSE41273_methylation_matrix.rds\n")

if (!file.exists(file.path(DATA_DIR, "GSE218676_metadata.csv"))) {
  suppressPackageStartupMessages(library(GEOquery))
  g <- getGEO("GSE218676", GSEMatrix = TRUE, getGPL = FALSE, destdir = DATA_DIR)
  pd <- dplyr::bind_rows(lapply(g, function(x) Biobase::pData(x)))
  rownames(pd) <- pd$geo_accession
  write.csv(pd, file.path(DATA_DIR, "GSE218676_metadata.csv"))
}

cat("\n=== STRING v12.0 human PPI ===\n")
fetch("https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz",
      file.path(DATA_DIR, "9606.protein.links.v12.0.txt.gz"), min_bytes = 1e7)

cat("\n=== TRRUST v2 human ===\n")
fetch("https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv",
      file.path(DATA_DIR, "trrust_rawdata.human.tsv"))

cat("\n=== Manifest ===\n")
files <- list.files(DATA_DIR, full.names = TRUE)
files <- files[!grepl("\\.part$|_series_matrix", files) & !dir.exists(files)]
manifest <- data.frame(file = basename(files),
                       bytes = file.size(files),
                       md5 = unname(tools::md5sum(files)),
                       fetched = format(Sys.time(), "%Y-%m-%d"))
write.csv(manifest, file.path(DATA_DIR, "MANIFEST.csv"), row.names = FALSE)
print(manifest[, c("file", "bytes")], row.names = FALSE)
cat("\nAll inputs ready.\n")
