# ============================================================
# Install every package the v5 pipeline loads.
# 1) install.packages() -> r2u/bspm apt binaries where they exist
# 2) BiocManager::install() for anything left (built into R_LIBS_USER,
#    which the notebook points at Google Drive so it persists)
# 3) GitHub-only packages (RRHO2)
# 4) verify every package loads; write versions to logs/
# ============================================================

options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 3600,
        Ncpus = max(1L, parallel::detectCores()))

user_lib <- Sys.getenv("R_LIBS_USER")
if (nzchar(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_lib, .libPaths()))
}
cat("R:", R.version.string, "\nLibraries:\n"); print(.libPaths())

cran <- c("BiocManager", "remotes", "tidyverse", "dplyr", "tidyr", "ggplot2",
          "ggrepel", "pheatmap", "reshape2", "scales", "igraph", "httr",
          "jsonlite", "ggvenn", "data.table", "viridisLite")
bioc <- c("DESeq2", "apeglm", "Biobase", "limma", "edgeR", "minfi",
          "IlluminaHumanMethylation450kanno.ilmn12.hg19", "EpiDISH", "GEOquery",
          "org.Hs.eg.db", "AnnotationDbi", "JASPAR2020", "TFBSTools", "motifmatchr",
          "BSgenome.Hsapiens.UCSC.hg38", "TxDb.Hsapiens.UCSC.hg38.knownGene",
          "GenomicFeatures", "GENIE3", "pandaR", "netZooR", "clusterProfiler",
          "decoupleR", "OmnipathR")
github <- c(RRHO2 = "Caleb-Huo/RRHO2")

have <- function(p) suppressWarnings(requireNamespace(p, quietly = TRUE))

# 1) binaries through bspm (fast). bspm handles both CRAN and Bioc names.
need <- c(cran, bioc)[!vapply(c(cran, bioc), have, logical(1))]
if (length(need)) {
  cat("\nBinary install attempt for", length(need), "packages...\n")
  try(install.packages(need), silent = FALSE)
}

# 2) whatever is still missing, from Bioconductor/CRAN source
need <- c(cran, bioc)[!vapply(c(cran, bioc), have, logical(1))]
if (length(need)) {
  cat("\nSource install for:", paste(need, collapse = ", "), "\n")
  if (!have("BiocManager")) install.packages("BiocManager")
  BiocManager::install(need, update = FALSE, ask = FALSE,
                       lib = if (nzchar(user_lib)) user_lib else .libPaths()[1])
}

# 3) GitHub
for (p in names(github)) if (!have(p)) {
  cat("\nGitHub install:", github[[p]], "\n")
  try(remotes::install_github(github[[p]], upgrade = "never", quiet = TRUE,
                              lib = if (nzchar(user_lib)) user_lib else .libPaths()[1]))
}

# 4) verify
all_pkgs <- c(cran, bioc, names(github))
status <- data.frame(package = all_pkgs,
                     loads = vapply(all_pkgs, have, logical(1)),
                     version = vapply(all_pkgs, function(p)
                       tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_),
                       character(1)))
print(status, row.names = FALSE)
dir.create(file.path(Sys.getenv("FXS_ROOT", "."), "logs"), showWarnings = FALSE, recursive = TRUE)
write.csv(status, file.path(Sys.getenv("FXS_ROOT", "."), "logs", "package_versions.csv"),
          row.names = FALSE)
cat("\nBioconductor:", as.character(BiocManager::version()), "\n")

bad <- status$package[!status$loads & status$package != "RRHO2"]
if (length(bad)) stop("These packages failed to install: ", paste(bad, collapse = ", "))
if (!status$loads[status$package == "RRHO2"])
  cat("RRHO2 missing — Step 16 falls back to the built-in RRHO implementation.\n")
cat("\nAll required packages load.\n")
