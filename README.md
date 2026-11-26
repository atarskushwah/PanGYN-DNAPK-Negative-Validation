# PanGYN-DNAPK-Negative-Validation

**Baseline Multi-Omic Profiles Fail to Predict Functional DNA-PK Dependency Across Cervical, Ovarian, and Endometrial Cancers: A Comparative Pan-Gynecologic Negative-Validation Study**

## Overview
This repository houses the reproducible R computational pipeline, harmonized summary tables, 600 DPI publication figures, and submission-ready manuscript files for a multi-lineage negative-validation audit of DNA-PK/c-NHEJ (`PRKDC`, `XRCC5`/Ku80, `XRCC6`/Ku70, `LIG4`, `NHEJ1`) dependency across **108 human gynecologic cancer cell lines**—Cervical (`CACX`, n = 18), Ovarian (`OVCX`, n = 56), and Endometrial (`ENDOCX`, n = 34)—benchmarked against **1,140 pan-cancer reference models** in DepMap Public 26Q1 and orthogonal **DEMETER2 v6 RNAi** screens.

### Key Findings
- **DepMap 26Q1 `PRKDC (5591)` Pseudogene Masking:** True CRISPR `PRKDC (5591)` is `NA`-masked in **1,134 / 1,208 (93.9%)** pan-cancer models and **106 / 108 (98.1%)** gynecologic models (`0/18` CACX, `2/56` OVCX, `0/34` ENDOCX) due to Broad Avana guide sequence homology with pseudogene *PRKDCP1*.
- **100% Holoenzyme & Orthogonal RNAi Coverage:** Obligatory DNA-PK holoenzyme partners `XRCC5` (Ku80), `XRCC6` (Ku70), `LIG4`, and `NHEJ1` have complete **108 / 108** CRISPR coverage alongside **52** orthogonal `DEMETER2` *PRKDC* RNAi profiles.
- **Zero Abundance-to-Dependency Coupling:** Baseline target mRNA expression shows no significant correlation with functional dependency in `CACX`, `OVCX`, `ENDOCX`, or `Pooled Pan-GYN` (all *P* > 0.06).
- **Universal Collapse to the Null Hypothesis:** Across cross-lineage pan-cancer transport, within-lineage nested cross-validation (19,210 leak-free mRNAs + WGS copy number + LoF somatic mutations), epithelial-only sensitivity cohorts (`N = 89`), and within-lineage `DEMETER2` *PRKDC* RNAi (`N = 52`), regularized models consistently fail to outperform the intercept-only null baseline ($\Delta\text{RMSE} \ge 0$).
- **REMARK & TRIPOD Clinical Halt:** Because predictive models collapse to zero-variance or non-generalizing null intercepts, downstream clinical projections into `TCGA-CESC` (n = 309), `TCGA-OV` (n = 379), `TCGA-UCEC` (n = 548), and GEO cohorts were formally halted.

## Repository Structure
```text
PanGYN-DNAPK-Negative-Validation/
├── .gitignore                                              # Excludes raw DepMap files & large .rds objects
├── README.md                                               # Project documentation & summary benchmarks
├── scripts/
│   ├── 01_pangyn_harmonization.R                           # DepMap 26Q1 + DEMETER2 harmonisation & Table 1
│   ├── 02_pangyn_comparative_cv.R                          # Leakage-free cross-lineage & nested CV (Table 2, S2A, S2B)
│   ├── 03_pangyn_codependency.R                            # Abundance vs. dependency & complex coupling (Table 3A, 3B, S3)
│   ├── 04_generate_pangyn_figures.R                        # Publication-grade 600 DPI Figures 1–4 (.tiff & .png)
│   └── 05_build_manuscript_docx.R                          # Programmatic Main + Supplementary .docx & README builder
├── results/
│   ├── tables/                                             # Verified CSV tables (Table 1–3, Table S1–S3)
│   └── figures/                                            # High-resolution 600 DPI figures
└── manuscript/
    ├── PanGYN_DNAPK_Negative_Validation_Main_v1.1.docx
    └── PanGYN_DNAPK_Negative_Validation_Supplementary_v1.1.docx
```

## Reproducing the Pipeline
1. Place DepMap Public 26Q1 (`Model.csv`, `CRISPRGeneEffect.csv`, `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`, `OmicsCNGeneWGS.csv`, `OmicsSomaticMutations.csv`) and `DEMETER2 v6` (`D2_combined_gene_dep_scores.csv`) inside `data/raw/`.
2. Execute `scripts/01_pangyn_harmonization.R` through `scripts/05_build_manuscript_docx.R` sequentially in R.
