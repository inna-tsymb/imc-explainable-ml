# install_packages.R
# ──────────────────────────────────────────────────────────────────────────────
# Fallback manual installation script.
# Preferred: use renv::restore() from renv.lock instead.
#
# Usage: source("install_packages.R")
# ──────────────────────────────────────────────────────────────────────────────

message("=== IMC XAI — Package Installation ===\n")

# ── 0. renv (environment manager) ────────────────────────────────────────────
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# ── 1. BiocManager ───────────────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(version = "3.18", ask = FALSE, update = FALSE)

# ── 2. Bioconductor packages ─────────────────────────────────────────────────
bioc_packages <- c(
  # Core single-cell data structures
  "SingleCellExperiment",
  "SpatialExperiment",
  "SummarizedExperiment",
  "S4Vectors",

  # Single-cell analysis
  "scater",       # quality control, dimensionality reduction
  "scuttle",      # aggregation, normalization utilities
  "scran",        # normalization, clustering support
  "batchelor",    # batch correction

  # IMC-specific
  "CATALYST",     # spillover correction, debarcoding
  "imcRtools",    # spatial graph, cell-cell interactions
  "cytomapper",   # image/mask loading, Shiny viewer

  # Visualization
  "dittoSeq",     # heatmaps, plots for single-cell data
  "ComplexHeatmap",

  # Parallelism
  "BiocParallel",
  "BiocNeighbors"
)

message("Installing Bioconductor packages...")
BiocManager::install(bioc_packages, ask = FALSE, update = FALSE)

# ── 3. CRAN packages ─────────────────────────────────────────────────────────
cran_packages <- c(
  # Machine learning
  "xgboost",    # gradient boosting classifier
  "caret",      # grouped CV, confusion matrix
  "shapviz",    # SHAP explanations for tree models

  # Tidyverse
  "tidyverse",  # ggplot2, dplyr, tidyr, stringr, purrr, forcats, readr
  "patchwork",  # combining ggplot figures

  # Visualization extras
  "viridis",
  "ggridges",
  "ggrepel",
  "RColorBrewer",
  "circlize",   # required by ComplexHeatmap
  "cowplot",

  # Utilities
  "doParallel",
  "umap",       # UMAP for image-level composition
  "rmarkdown",
  "knitr"
)

message("Installing CRAN packages...")
install.packages(
  cran_packages,
  repos = "https://cloud.r-project.org",
  dependencies = TRUE
)

# ── 4. Verify installation ───────────────────────────────────────────────────
all_packages <- c(bioc_packages, cran_packages)
all_packages <- unique(all_packages[!all_packages %in% c("tidyverse")])

# tidyverse is a meta-package — check its components individually
tidyverse_core <- c("ggplot2", "dplyr", "tidyr", "stringr", "purrr", "forcats")
check_packages <- c(all_packages, tidyverse_core)

missing <- check_packages[
  !sapply(check_packages, requireNamespace, quietly = TRUE)
]

if (length(missing) == 0) {
  message("\n✅  All packages installed successfully.")
} else {
  message("\n⚠️  The following packages could not be installed:")
  message(paste(" -", missing, collapse = "\n"))
  message('\nTry: BiocManager::install(c("', paste(missing, collapse = '", "'), '"))')
}

message("\nR version: ", R.version.string)
message("Bioconductor version: ", BiocManager::version())
