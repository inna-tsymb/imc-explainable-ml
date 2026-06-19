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
| `03a_uncertainty_analysis`, `03b_uncertainty_comparison` | **Uncertainty** | Which cells should we not trust? |
| `04_spatial_xai` | **Spatial context** | Do neighbours influence the label? |
| `05_predict_unlabeled` | **Unlabeled images** | Obtain predictions for images with no manual labels |
| `05_visual_model_validation` | **Image-level** | Do model predictions visually make sense? |
| `05_dataset_xgboost_composition_analysis` | **Dataset-level** | Do predicted cell types make sense composition-wise? |
| `05_anomalies_handoff` | **Image-level** | PCA analysis of samples to identify possible outliers |

---

## Dataset

- **Cohort:** Colorectal cancer (CRC) tissue — University Saint-Luc (USL), Brussels
- **Panel:** Panel_2_10 · 44 protein markers
- **Labeled images:** 27 IMC acquisitions used for model training and evaluation (1 blank image removed after QC)
- **Unseen images:** 44 additional IMC acquisitions without manual labels — cell type labels obtained exclusively via the trained XGBoost classifier (`05_predict_unlabeled.Rmd`)
- **Total dataset:** 71 images · ~250,000 single cells
- **Cell types (20):** B, BnT, CD4, CD8, Treg, NK, Igg, DC, Mac, MacCD204, MacCD209, mregDC, Neutrophil, PMN_MDSC, MDSC, Vasculature, vCAF, Tumor, SMA, Fibro, undefined

> Raw data are **not included** in this repository.
> Place steinbock output under `data/Panel_2_10/` (see [Data Setup](#data-setup)).

---

## Repository Structure

```
.
├── notebooks/
│   ├── 01_read_data.Rmd                             # Preprocessing, spillover, gating, XGBoost CV + final model
│   ├── 02_shap_analysis.Rmd                         # SHAP values, feature importance, beeswarm, interaction plots
│   ├── 03a_uncertainty_analysis.Rmd                 # Softmax entropy, per-class uncertainty on labeled test set
│   ├── 03b_uncertainty_comparison.Rmd               # Uncertainty comparison: labeled test set vs 44 unseen images
│   ├── 04_spatial_xai.Rmd                           # Spatial graph, neighbour composition, cell-cell interactions
│   ├── 05_predict_unlabeled.Rmd                     # XGBoost predictions on 44 unseen images
│   ├── 05_visual_model_validation.Rmd               # Visual validation of Tumor and mregDC predictions
│   ├── 05_dataset_xgboost_composition_analysis.Rmd  # CLR-transformed composition heatmap across all 71 images
│   └── 05_anomalies_handoff.Rmd                     # PCA outlier detection at sample level
├── results/
│   ├── figures/                                     # PNG/PDF outputs from all notebooks
│   └── tables/                                      # .csv outputs from all notebooks
├── install_packages.R                               # R package installation script
├── LICENSE
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
01_read_data.Rmd (R / Bioconductor)           [27 labeled images]
  ├── read_steinbock() → SpatialExperiment
  ├── Spillover correction (CATALYST)
  ├── Manual gating labels (cytomapper)
  ├── Label consolidation + QC
  ├── Stratified train/test split (image-level)
  ├── XGBoost grid search (5-fold grouped CV)
  └── Final model → classifier_xgboost.rds
          ↓                                   ↓
  XAI on labeled data             05_predict_unlabeled.Rmd
  ┌─────────────────┐             [44 unseen images → ~250k cells]
  │ 02  SHAP        │                     ↓
  │ 03a Uncertainty │             05_visual_model_validation.Rmd
  │ 04  Spatial XAI │             05_dataset_xgboost_composition_analysis.Rmd
  └─────────────────┘             05_anomalies_handoff.Rmd
          ↓                                   ↓
  03b_uncertainty_comparison.Rmd  [labeled test set vs unseen images]
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

Top markers by XGBoost gain: `CD3`, `CD204`, `CD209`, `CD20`, `Lamp3`, `CD31`, `FOXP3`, `CD66b`, `CD146`.

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

Or run the provided script: `source("install_packages.R")`

---

## Key Outputs

### Model validation and predictions on unseen data (`05_*.Rmd`)

| File | Description |
|------|-------------|
| `01_Validation_1A_Predicted_TME.pdf` | Spatial overlay of XGBoost-predicted cell types for 2 labeled training images and 2 unseen images, coloured by cell type — visual sanity check that predictions match known tissue architecture |
| `02_Validation_1B_GroundTruth_TME.pdf` | Same spatial overlay using manual gating labels for the 2 training images — used as ground truth reference for comparison with predicted TME |
| `03_Validation_1C_HighContrast_Tumor.pdf` | Tumor cells only, displayed in high contrast on the same images — validates that predicted tumor shapes align with expected epithelial structures |
| `04_Validation_2A_Predicted_Tumor_Overlay.pdf` | Predicted tumor cells overlaid on a high-contrast composite of panCK and Ecad marker channels — confirms marker-to-cell-type correspondence for the classifier |
| `05_Validation_2B_GroundTruth_Tumor_Overlay.pdf` | Same panCK/Ecad overlay using original manual labels — paired comparison to `04` |
| `06_Validation_3_Predicted_mregDC.pdf` | Predicted mregDC cells overlaid on panCK and Ecad channels — validates that mregDC predictions do not co-localise with epithelial markers, as expected biologically |
| `07a_RAW_Proportions_Per_Image.pdf` | Stacked bar chart of raw predicted cell type proportions per image across all 71 images — reveals inter-image variability in TME composition |
| `07b_QC_Dead_Images_Marker_Heatmap.pdf` | Marker expression heatmap restricted to images with a high proportion of undefined predictions — used to identify low-quality acquisitions or staining failures |
| `08_CLEAN_CLR_Phenotype_Heatmap.pdf` | Centred log-ratio (CLR) transformed cell type proportions per image, displayed as a clustered heatmap — reveals sample-level groupings and TME subtypes across the full 71-image dataset |
| `09_PCA_Anomaly_Map.pdf` | PCA of CLR-transformed image compositions with Mahalanobis-distance outlier scores — identifies atypical images for downstream biological investigation |

### SHAP explainability (`02_shap_analysis.Rmd`)

| File | Description |
|------|-------------|
| `shap_global_importance.png` | Horizontal bar chart of mean absolute SHAP values per marker, averaged across all cell types and test cells — global ranking of marker contribution to classifier decisions |
| `shap_beeswarm_global.png` | SHAP beeswarm plot across all test cells: each dot is one cell, x-axis is SHAP value, colour encodes marker expression level — shows how high vs low expression of each marker pushes predictions |
| `shap_class_heatmap.png` | Heatmap of mean absolute SHAP values per marker per predicted cell type — reveals which markers are most diagnostic for each cell type |
| `shap_interaction_CD3_CD8a.png` | SHAP dependence plot showing how the contribution of CD3 varies with CD8a expression level — confirms that the model distinguishes CD4 and CD8 T cells through marker interaction |
| `shap_blur_SMA_CD146.png` | Dependence plot illustrating how co-expression of CD146 and SMA shifts predictions between vCAF and SMA cell types — highlights the decision boundary between two stromal populations |

### Uncertainty quantification (`03a_uncertainty_analysis.Rmd`, `03b_uncertainty_comparison.Rmd`)

| File | Description |
|------|-------------|
| `uncertainty_distributions.png` | Overlapping histograms of max softmax probability and normalised prediction entropy for the labeled test set, split by correct vs incorrect predictions |
| `uncertainty_ridgeplot.png` | Ridge plot of prediction entropy distributions per cell type — ranks cell types by median uncertainty and reveals which populations are hardest to classify |
| `uncertainty_confusion_pairs.png` | Bar chart of the top-20 most frequent misclassification pairs on the test set — shows which cell types are most often confused with one another |
| `uncertainty_entropy_margin_scatter.png` | Scatter plot of entropy vs probability margin (P1 − P2), coloured by prediction correctness — identifies cells that are both uncertain and at high risk of misclassification |
| `uncertainty_comparison_*.png` | Comparison panels between labeled test set and 44 unseen images: entropy distributions, per-class uncertainty rates, and high-uncertainty cell fractions |

### Spatial XAI (`04_spatial_xai.Rmd`)

| File | Description |
|------|-------------|
| `s01_interaction_heatmap` | Heatmap of average z-scores for all cell-type pair interactions across the 27 labeled images, with dendrogram clustering of cell types by spatial co-occurrence patterns |
| `s02_cd8_tumor_proximity` | Per-image histogram of CD8–Tumour spatial proximity z-scores — quantifies how consistently CD8 T cells are enriched in the neighbourhood of tumour cells across samples |
| `s03_tumor_neighbour_fraction` | Boxplots showing the fraction of direct 15 µm neighbours that are Tumour cells, broken down by the cell type of the focal cell — reveals which cell types preferentially border tumour regions |
| `s04_augmented_feature_importance` | Bar chart of the top 30 features in an augmented XGBoost model trained on both marker expression and spatial neighbourhood composition — quantifies the additional predictive value of spatial context |

---

## Contributors

**Uliana Krektun**
- `02_shap_analysis.Rmd` — SHAP feature importance, beeswarm, interaction plots (CD3×CD8a, CD146×SMA blur), SHAP class heatmap
- Outputs: `shap_global_importance.png`, `shap_beeswarm_global.png`, `shap_class_heatmap.png`, `shap_interaction_CD3_CD8a.png`, `shap_blur_SMA_CD146.png`

**Inna Kucherova**
- `03a_uncertainty_analysis.Rmd`, `03b_uncertainty_comparison.Rmd` — softmax entropy, per-class uncertainty, confusion pairs, prediction confidence distributions, labeled vs unseen comparison
- Repository setup: README, `.gitignore`, `install_packages.R`, repository structure and description

**Anna Gurina**
- `04_spatial_xai.Rmd` — spatial neighbour graph, cell-cell interaction heatmap, CD8–Tumour proximity, augmented model with spatial features
- Outputs: `s01_interaction_heatmap`, `s02_cd8_tumor_proximity`, `s03_tumor_neighbour_fraction`, `s04_augmented_feature_importance`

**Oleksandra Tsepilova**
- `05_predict_unlabeled.Rmd`, `05_visual_model_validation.Rmd`, `05_dataset_xgboost_composition_analysis.Rmd`, `05_anomalies_handoff.Rmd` — predictions on 44 unseen images, visual validation (Tumor, mregDC), CLR-transformed composition heatmap, PCA outlier detection
- Outputs: `01_Validation_1A` through `06_Validation_3`, `07a/b`, `08_CLEAN_CLR_Phenotype_Heatmap.pdf`, `09_PCA_Anomaly_Map.pdf`

**All team members**
- `01_read_data.Rmd` — preprocessing, spillover correction, manual gating, XGBoost cross-validation and final model training
- `01b_save_shared_objects.Rmd` — shared RDS objects across all notebooks
- Project concept development: definition of research questions, XAI workstream design, and division of analytical tasks
- Repository documentation: README, LICENSE, `.gitignore`, `install_packages.R`
- Presentation design, structure, and delivery
- Joint interpretation of results and conclusions

---

## Responsible Use of AI

Claude (Anthropic) and Gemini (Google) were used during this project for the following purposes:

- **Debugging** — identifying and resolving errors in R code across analysis notebooks
- **Concept explanation** — clarifying IMC-specific concepts, Bioconductor data structures, and XAI methodology (SHAP, entropy, spatial graph analysis)
- **Brainstorming** — discussing analytical approaches, workstream division, and interpretation of results
- **Code standardisation** — ensuring consistent style, explicit namespace prefixes (`dplyr::`, `base::`), and translation of code comments to English across all notebooks
- **Writing style** — editing prose in the README, repository description, and presentation materials to match a formal academic register

All analytical decisions, biological interpretations, and final code were reviewed and validated by the team members.

---

## Citation

```
@misc{imc_xai_crc_2026,
  title  = {Explainable ML for Cellular Phenotyping in Imaging Mass Cytometry},
  author = {Uliana Krektun, Inna Kucherova, Anna Gurina, Oleksandra Tsepilova},
  year   = {2026},
  url    = {https://github.com/inna-tsymb/imc-explainable-ml}
}
```

---

## License

MIT — see [LICENSE](LICENSE).

---
