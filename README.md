# Explainable ML for Cellular Phenotyping in Imaging Mass Cytometry

[![R 4.3+](https://img.shields.io/badge/R-4.3%2B-276DC3?logo=r)](https://www.r-project.org/)
[![Bioconductor 3.18](https://img.shields.io/badge/Bioconductor-3.18-brightgreen)](https://bioconductor.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

This repository contains the full analysis pipeline for **cellular phenotyping** in
**Imaging Mass Cytometry (IMC)** data from colorectal cancer (CRC) tissue sections,
with a focus on **Explainable AI (XAI)**.

A multiclass **XGBoost** classifier is trained on manually gated single-cell IMC
data and then interrogated through four complementary XAI lenses:

| Notebook | Focus | Key question |
|----------|-------|--------------|
| `02_shap_analysis` | **SHAP** | Why did the model assign this cell type? |
| `03_uncertainty_analysis` | **Uncertainty** | Which cells should we not trust? |
| `04_spatial_xai` | **Spatial context** | Do neighbours influence the label? |
| `05_image_composition_xai` | **Image-level** | Which images are typical vs. anomalous? |

---

## Dataset

- **Cohort:** Colorectal cancer (CRC) tissue — University Saint-Luc (USL), Brussels
- **Panel:** Panel_2_10 (~40 protein markers)
- **Images:** 27 IMC acquisitions (1 blank image removed after QC)
- **Cell types (20):** B, BnT, CD4, CD8, Treg, NK, Igg, DC, Mac, MacCD204,
  MacCD209, mregDC, Neutrophil, PMN_MDSC, MDSC, Vasculature, vCAF, Tumor, SMA,
  Fibro, undefined

> Raw data are **not included** in this repository.  
> Place steinbock output under `data/Panel_2_10/` (see [Data Setup](#data-setup)).

---

## Repository Structure

```
.
├── notebooks/
│   ├── 01_read_data.Rmd              # Preprocessing, spillover, gating, XGBoost CV + final model
│   ├── 01b_save_shared_objects.Rmd   # Saves shared RDS objects for notebooks 02–05
│   ├── 02_shap_analysis.Rmd          # SHAP values, beeswarm, waterfall, interaction plots
│   ├── 03_uncertainty_analysis.Rmd   # Softmax entropy, per-class uncertainty, confusion pairs
│   ├── 04_spatial_xai.Rmd            # Spatial graph, neighbour composition, cell interactions
│   └── 05_image_composition_xai.Rmd  # PCA/UMAP of images, TME subtypes, outlier detection
├── data/
│   └── Panel_2_10/                   # steinbock output (not tracked by git)
│       ├── img/
│       ├── masks/
│       ├── intensities/
│       ├── regionprops/
│       ├── neighbors/
│       ├── compensation/
│       └── gates/
├── results/
│   ├── figures/                      # PNG outputs from all notebooks
│   └── tables/                       # CSV summaries
├── renv.lock                         # Pinned R environment (renv)
├── .Rprofile                         # Auto-activates renv
└── README.md
```

---

## Pipeline Overview

```
steinbock CLI (Python)
  ├── preprocess imc images --hpf 50
  ├── segment deepcell --minmax
  ├── measure intensities / regionprops
  └── measure neighbors --type centroids --dmax 15
          ↓
01_read_data.Rmd (R / Bioconductor)
  ├── read_steinbock() → SpatialExperiment
  ├── Spillover correction (CATALYST)
  ├── Manual gating labels (cytomapper)
  ├── Label consolidation + QC
  ├── Stratified train/test split (image-level)
  ├── XGBoost grid search (5-fold grouped CV)
  └── Final model → classifier_xgboost.rds
          ↓
01b_save_shared_objects.Rmd
  └── classifier_inputs.rds  (X, test_prob, label_mapping, …)
          ↓
┌─────────────────────────────────────────────────┐
│  02  SHAP          │  03  Uncertainty           │
│  04  Spatial XAI   │  05  Image composition     │
└─────────────────────────────────────────────────┘
```

---

## Classifier Summary

| Parameter | Value |
|-----------|-------|
| Model | XGBoost `multi:softprob` |
| `max_depth` | 4 |
| `eta` | 0.1 |
| `subsample` | 0.8 |
| `colsample_bytree` | 0.8 |
| `min_child_weight` | 5 |
| `nrounds` | 108 (early stopping) |
| CV strategy | 5-fold grouped (by `sample_id`) |
| Overall accuracy (test) | **95.1%** |
| Cohen's κ (test) | **0.896** |

Top markers by XGBoost gain: `CD3`, `CD204`, `CD209`, `CD20`, `Lamp3`, `CD31`,
`FOXP3`, `CD66b`, `CD146`.

---

## Data Setup

1. Run steinbock on your raw IMC data:
```bash
steinbock preprocess imc images --hpf 50
steinbock segment deepcell --minmax
steinbock measure intensities
steinbock measure regionprops
steinbock measure neighbors --type centroids --dmax 15
```

2. Place the output in `data/Panel_2_10/`

3. Place gating `.rds` files from cytomapper in `data/Panel_2_10/gates/`

4. Update `base_path` in `01_read_data.Rmd` if your directory differs.

---

## Environment Setup

### Option A — renv (recommended)

```r
# In R, from the project root:
install.packages("renv")
renv::restore()   # installs all packages from renv.lock
```

### Option B — manual install

```r
# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "SingleCellExperiment", "SpatialExperiment", "SummarizedExperiment",
  "scater", "scuttle", "scran", "CATALYST", "imcRtools",
  "cytomapper", "dittoSeq", "BiocNeighbors", "BiocParallel", "batchelor"
))

# CRAN packages
install.packages(c(
  "xgboost", "caret", "shapviz",
  "ggplot2", "dplyr", "tidyr", "tidyverse", "forcats",
  "patchwork", "viridis", "ggridges", "ggrepel", "RColorBrewer",
  "ComplexHeatmap", "circlize",
  "stringr", "purrr", "doParallel", "cowplot", "umap"
))
```

See `renv.lock` for exact pinned versions.

---

## Running the Analysis

```r
# 1. Render the main preprocessing + classifier notebook
rmarkdown::render("notebooks/01_read_data.Rmd")

# 2. Save shared objects (add this block to 01_read_data.Rmd or run manually)
rmarkdown::render("notebooks/01b_save_shared_objects.Rmd")

# 3. Run XAI notebooks independently (in any order)
rmarkdown::render("notebooks/02_shap_analysis.Rmd")
rmarkdown::render("notebooks/03_uncertainty_analysis.Rmd")
rmarkdown::render("notebooks/04_spatial_xai.Rmd")
rmarkdown::render("notebooks/05_image_composition_xai.Rmd")
```

All figures are saved to `results/figures/` and tables to `results/tables/`.

---

## Key Outputs

| File | Notebook | Description |
|------|----------|-------------|
| `shap_global_importance.png` | 02 | Global mean \|SHAP\| bar chart |
| `shap_beeswarm_global.png` | 02 | SHAP beeswarm — all test cells |
| `shap_beeswarm_per_class.png` | 02 | Per-class beeswarm (CD8, CD4, Tumor, …) |
| `shap_class_heatmap.png` | 02 | Marker × cell-type SHAP heatmap |
| `shap_interaction_CD3_CD8a.png` | 02 | CD3 × CD8a dependence plot |
| `uncertainty_distributions.png` | 03 | max_prob and entropy distributions |
| `uncertainty_ridgeplot.png` | 03 | Entropy by cell type (ridge plot) |
| `uncertainty_confusion_pairs.png` | 03 | Top-20 confusion pairs |
| `spatial_neighbor_composition.png` | 04 | Neighbour composition heatmap |
| `spatial_cell_interactions.png` | 04 | Cell-cell enrichment (permutation test) |
| `spatial_uncertainty_map.png` | 04 | Spatial map of prediction entropy |
| `image_composition_heatmap.png` | 05 | Image × cell-type heatmap |
| `image_composition_pca.png` | 05 | PCA of image compositions |
| `image_outlier_detection.png` | 05 | Anomalous images (Mahalanobis distance) |

---

## Citation

If you use this code, please cite:

```
@misc{imc_xai_crc_2024,
  title  = {Explainable ML for Cellular Phenotyping in Imaging Mass Cytometry},
  author = {<Your Name>},
  year   = {2024},
  url    = {https://github.com/<your-username>/<repo-name>}
}
```

---

## License

MIT — see [LICENSE](LICENSE).
