# PanGYN-DNAPK-Negative-Validation

**Baseline Multi-Omic Profiles Fail to Predict Functional DNA-PK Dependency Across Cervical, Ovarian, and Endometrial Cancers: A Comparative Pan-Gynecologic Negative-Validation Study**

- **Author:** Tanmay Singh Naruka (Department of Biotechnology, Amity University, Noida, Uttar Pradesh, India)
- **Repository:** https://github.com/atarskushwah/PanGYN-DNAPK-Negative-Validation

## Overview
This repository houses the reproducible R computational pipeline, harmonized summary tables, and 600 DPI publication figures for a multi-lineage negative-validation audit of DNA-PK/c-NHEJ (`PRKDC`, `XRCC5`/Ku80, `XRCC6`/Ku70, `LIG4`, `NHEJ1`) dependency across **108 human gynecologic cancer cell lines**—Cervical (`CACX`, n = 18), Ovarian (`OVCX`, n = 56), and Endometrial (`ENDOCX`, n = 34)—benchmarked against **1,140 pan-cancer reference models** in DepMap Public 26Q1 and orthogonal **DEMETER2 v6 RNAi** screens.

### Key Findings
- **DepMap 26Q1 `PRKDC (5591)` Pseudogene Masking:** True CRISPR `PRKDC (5591)` is `NA`-masked in **1,134 / 1,208 (93.9%)** pan-cancer models and **106 / 108 (98.1%)** gynecologic models (`0/18` CACX, `2/56` OVCX, `0/34` ENDOCX) due to Broad Avana guide sequence homology with pseudogene *PRKDCP1*.
- **100% Holoenzyme & Orthogonal RNAi Coverage:** Obligatory DNA-PK holoenzyme partners `XRCC5` (Ku80), `XRCC6` (Ku70), `LIG4`, and `NHEJ1` have complete **108 / 108** CRISPR coverage alongside **52** orthogonal `DEMETER2` *PRKDC* RNAi profiles.
- **Zero Abundance-to-Dependency Coupling:** Baseline target mRNA expression shows no significant correlation with functional dependency in `CACX`, `OVCX`, `ENDOCX`, or `Pooled Pan-GYN` (all *P* > 0.06).
- **Universal Collapse to the Null Hypothesis:** Across cross-lineage pan-cancer transport, within-lineage nested cross-validation (19,210 leak-free mRNAs + WGS copy number + LoF somatic mutations), epithelial-only sensitivity cohorts (`N = 89`), and `DEMETER2` *PRKDC* RNAi (`N = 52`), regularized models consistently fail to outperform the intercept-only null baseline (ΔRMSE ≥ 0).
- **REMARK & TRIPOD Clinical Halt:** Because predictive models collapse to zero-variance or non-generalizing null intercepts, downstream clinical projections into `TCGA-CESC` (n = 309), `TCGA-OV` (n = 379), `TCGA-UCEC` (n = 548), and GEO cohorts were formally halted.

## Repository Contents

### Reproducible R Scripts
- `01_pangyn_harmonization.R` — DepMap Public 26Q1 & DEMETER2 v6 harmonization, identifier audit, and Table 1 generation
- `02_pangyn_comparative_cv.R` — Leakage-free cross-lineage and repeated nested cross-validation benchmark (Tables 2, S2A, S2B)
- `03_pangyn_codependency.R` — Baseline mRNA abundance vs. dependency and lineage-adjusted complex co-dependency (Tables 3A, 3B, S3)
- `04_generate_pangyn_figures.R` — Publication-grade 600 DPI Figures 1–4
- `05_build_manuscript_docx.R` — Programmatic manuscript document builder

### Harmonized Summary Tables (`.csv`)
- `Table1_PanGYN_Cohort_Architecture.csv` — Discovery cohort architecture & assay completeness (`N = 108`)
- `Table2_PanGYN_Nested_CV_Benchmark.csv` — Comparative cross-lineage and within-lineage nested CV benchmark
- `Table3A_Abundance_vs_Dependency_PanGYN.csv` — Univariate mRNA abundance vs. perturbation dependency correlations
- `Table3B_DNAPK_Holoenzyme_CoDependency.csv` — Lineage-adjusted co-dependency with `XRCC5` (Ku80)
- `TableS1_PanGYN_108_CellLines_Manifest.csv` — Full manifest of all 108 gynecologic cell line models
- `TableS2A_Epithelial_Sensitivity_CV.csv` — Epithelial-only sensitivity nested CV benchmark (`N = 89`)
- `TableS2B_DEMETER2_PRKDC_RNAi_CV.csv` — Orthogonal DEMETER2 *PRKDC* RNAi benchmark (`N = 52`)
- `TableS3_DEMETER2_RNAi_vs_CRISPR_Concordance.csv` — CRISPR vs. DEMETER2 RNAi perturbation concordance

### Publication Figures (600 DPI `.png`)
- `Figure1_PanGYN_Cohort_and_Assay_Audit.png`
- `Figure2_PanGYN_Dependency_and_Abundance.png`
- `Figure3_PanGYN_NestedCV_Null_Collapse.png`
- `Figure4_PanGYN_Holoenzyme_CoDependency.png`

## Reproducing the Pipeline
1. Download DepMap Public 26Q1 (`Model.csv`, `CRISPRGeneEffect.csv`, `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`, `OmicsCNGeneWGS.csv`, `OmicsSomaticMutations.csv`) and DEMETER2 v6 (`D2_combined_gene_dep_scores.csv`) from the Broad Institute DepMap portal into `data/raw/`.
2. Run `01_pangyn_harmonization.R` through `05_build_manuscript_docx.R` sequentially in R.

## Citation
If you use these workflows or benchmark tables, please cite `https://github.com/atarskushwah/PanGYN-DNAPK-Negative-Validation` alongside the Broad Institute DepMap Public 26Q1 and DEMETER2 v6 releases.
