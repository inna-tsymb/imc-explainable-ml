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
| `05_predict_unlabeled` | **Unlabeled Images** | Obtain prediction for images, that have no manual labels (therefore not present in train or test sets) |
| `05_visual_model_validation` | **Image-Level** | Check if model predictions visually make sense. Closer attention is paid to tumor and mregDC cell types. |
| `05_dataset_xgboost_composition_analysis` | **Dataset-Level** | Do predicted cell types make sense composition-wise for this dataset? |
`05_anomalies_handoff` | **Image-Level** | PCA analysis of samples to understand possible outliers and pass them to other checks.|

---

## Dataset

- **Cohort:** Colorectal cancer (CRC) tissue — University Saint-Luc (USL), Brussels
- **Panel:** Panel_2_10 (44 protein markers)
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
│   ├── 02_shap_analysis.Rmd          # SHAP values, feature importance bar chart, beeswarm, interaction plots
│   ├── 03_uncertainty_analysis.Rmd   # Softmax entropy, per-class uncertainty, confusion pairs
│   ├── 04_spatial_xai.Rmd            # Spatial graph, neighbour composition, cell interactions
│   └── 05_<...>.Rmd  # Unseen data analysis and validation, PCA/UMAP of images, TME subtypes, outlier detection
├── results/
│   └── figures/                       # PNG/PDF and other plots from script
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

---

## Key Outputs

| File | Description |
|------|-------------|
| `01_Validation_1A_Predicted_TME.pdf` | Visualization of predicted cell types for 2 training and 2 unseen images |
| `02_Validation_1B_GroundTruth_TME.pdf` | Visualization of ground truth labels for the same 2 training images |
| `03_Validation_1C_HighContrast_Tumor.pdf` | Visualization of gsame samples but only tumor cells, for validation of tumor shapes |
| `04_Validation_2A_Predicted_Tumor_Overlay.pdf` | Predicted tumor cells overlayed on top of high-contrast display of panCK and Ecad markers, for marker-cell type validation |
| `05_Validation_2B_GroundTruth_Tumor_Overlay.pdf` | Same visuals applied to original labels |
| `06_Validation_3_Predicted_mregDC.pdf` | Predicted mregCD cells overlayed on top of high-contrast display of panCK and Ecad markers, for marker-cell type validation |
| `07_RAW_Proportions_Per_Image.pdf` | Proportions of predicted cell types per image |
| `08_CLEAN_CLR_Phenotype_Heatmap.pdf` | CLR-transforme predicted cell types per image, displayed on heatmap to see presence of sample-level connections  |
| `09_PCA_Anomaly_Map.pdf` | Visualization of outliers detection (sample-level) for further analysis  |
| `shap_global_importance.png` | Bar chart that shows the importance of markers for each cell type |
| `shap_beeswarm_global.png` | SHAP beeswarm — summary of how a XGBoost makes decisions across all test cells |
| `shap_class_heatmap.png` | Marker × cell-type SHAP heatmap |
| `shap_interaction_CD3_CD8a.png` | CD3 × CD8a dependence plot |
| `shap_myeloid_comparison.png` | Comparison of the contribution of markers to distinguishing between MacCD204 and MacCD209 |
| `shap_treg_foxp3_threshold.png` | Representation of how the model determines the structural identity of Tregs using a combination of the FOXP3 and CD3 markers |
| `s01_interaction_heatmap` | Interaction heatmap: each cell shows the average z-score for that cell-type pair across all 13 images. The dendrogram clusters cell types with similar spatial patterns.
| `s02_cd8_tumor_proximity` | CD8–Tumour proximity: A histogram where each bar = one image, x-axis = z-score for that image.  |
| `s03_tumor_neighbour_fraction` | Tumour neighbour fraction, answering what fraction of each cell direct 15 µm neighbours are Tumour cells? Shown as boxplots per cell type. |
| `s04_augmented_feature_importance` |Augmented feature importance: Top 30 features in augmented model (markers vs spatial)|


#TODO check what relevant after here

| `uncertainty_distributions.png` | max_prob and entropy distributions |
| `uncertainty_ridgeplot.png` | Entropy by cell type (ridge plot) |
| `uncertainty_confusion_pairs.png` |  Top-20 confusion pairs |
| `spatial_neighbor_composition.png` | Neighbour composition heatmap |
| `spatial_cell_interactions.png` | Cell-cell enrichment (permutation test) |
| `spatial_uncertainty_map.png` | Spatial map of prediction entropy |
| `image_composition_heatmap.png` | Image × cell-type heatmap |
| `image_composition_pca.png` | PCA of image compositions |
| `image_outlier_detection.png` | Anomalous images (Mahalanobis distance) |

---

## Citation

If you use this code, please cite:

```
@misc{imc_xai_crc_2026,
  title  = {Explainable ML for Cellular Phenotyping in Imaging Mass Cytometry},
  author = {Uliana Krektun, Inna Kucherova, Anna Gurina, Oleksandra Tsepilova},
  year   = {2026},
  url    = {https://github.com/inna-tsymb/imc-explainable-ml/tree/main}
}
```

---

## License

MIT — see [LICENSE](LICENSE).
