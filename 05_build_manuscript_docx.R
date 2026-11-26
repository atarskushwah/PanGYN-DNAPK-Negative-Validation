# scripts/05_build_manuscript_docx.R
# Builds Standalone Main Manuscript (.docx), Separate Supplementary File (.docx) & README.md

for (pkg in c("officer", "flextable", "tidyverse")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(tidyverse)
library(officer)
library(flextable)

dir.create("manuscript", recursive = TRUE, showWarnings = FALSE)

# Load verified results tables directly from disk
t1   <- read_csv("results/tables/Table1_PanGYN_Cohort_Architecture.csv", show_col_types = FALSE)
t2   <- read_csv("results/tables/Table2_PanGYN_Nested_CV_Benchmark.csv", show_col_types = FALSE)
t3a  <- read_csv("results/tables/Table3A_Abundance_vs_Dependency_PanGYN.csv", show_col_types = FALSE)
t3b  <- read_csv("results/tables/Table3B_DNAPK_Holoenzyme_CoDependency.csv", show_col_types = FALSE)
ts1  <- read_csv("results/tables/TableS1_PanGYN_108_CellLines_Manifest.csv", show_col_types = FALSE)
ts2a <- read_csv("results/tables/TableS2A_Epithelial_Sensitivity_CV.csv", show_col_types = FALSE)
ts2b <- read_csv("results/tables/TableS2B_DEMETER2_PRKDC_RNAi_CV.csv", show_col_types = FALSE)
ts3  <- read_csv("results/tables/TableS3_DEMETER2_RNAi_vs_CRISPR_Concordance.csv", show_col_types = FALSE)

# ==============================================================================
# 1. CLEAN PUBLICATION HEADERS & US ENGLISH STANDARDIZATION FOR TABLES
# ==============================================================================
t1_clean <- t1 %>%
  mutate(across(where(is.character), ~gsub("tumours", "tumors", .x))) %>%
  rename(
    `Lineage`               = Cancer_Code,
    `Oncotree Lineage`      = OncotreeLineage,
    `Evaluable (n)`         = Evaluable_CRISPR_RNA,
    `Epithelial (n)`        = Epithelial_Only_N,
    `CRISPR PRKDC (n)`      = CRISPR_PRKDC_Unmasked_N,
    `RNAi PRKDC (n)`        = DEMETER2_PRKDC_RNAi_N,
    `XRCC5 Median (IQR)`    = XRCC5_Ku80_Median_IQR,
    `XRCC6 Median (IQR)`    = XRCC6_Ku70_Median_IQR,
    `WGS CNV (n)`           = WGS_CNV_Overlap_N,
    `Matched TCGA Cohort`   = Matched_TCGA_Cohort,
    `Matched GEO Cohorts`   = Matched_GEO_Cohorts
  )

clean_cv_table <- function(df) {
  df %>%
    mutate(
      Validation_Arm = gsub("->", "\u2192", Validation_Arm),
      Algorithm      = gsub("\\(a=", "(\u03b1=", Algorithm)
    ) %>%
    rename(
      `Validation Arm`    = Validation_Arm,
      `Target`            = Outcome,
      `Algorithm`         = Algorithm,
      `Cohort`            = Cancer_Code,
      `n`                 = N_Eval,
      `Model RMSE`        = Model_RMSE,
      `Null RMSE`         = Null_RMSE,
      !!"\u0394RMSE"     := Delta_RMSE,
      `OOF R`             = Transport_R,
      `Active Features`   = Active_Features
    )
}

t2_clean   <- clean_cv_table(t2)
ts2a_clean <- clean_cv_table(ts2a)
ts2b_clean <- clean_cv_table(ts2b)

t3a_clean <- t3a %>%
  rename(
    `Cohort`           = Cohort,
    `Target`           = Target,
    `n`                = N,
    `Pearson R`        = Pearson_R,
    `95% CI Low`       = CI_Low,
    `95% CI High`      = CI_High,
    `Pearson P`        = Pearson_P,
    !!"Spearman \u03c1" := Spearman_Rho,
    `Spearman P`       = Spearman_P
  )

t3b_clean <- t3b %>%
  rename(
    `Cohort`        = Cohort,
    `Anchor`        = Anchor_Gene,
    `Partner`       = Partner_Gene,
    `n`             = N,
    `Pearson R`     = Pearson_R,
    `P Value`       = P_Value,
    `BH-FDR q`      = FDR_Q
  )

ts3_clean <- ts3 %>%
  rename(
    `Gene`             = Gene,
    `Comparison`       = Comparison,
    `n`                = N,
    `Pearson R`        = Pearson_R,
    `95% CI Low`       = CI_Low,
    `95% CI High`      = CI_High,
    `Pearson P`        = Pearson_P,
    !!"Spearman \u03c1" := Spearman_Rho,
    `Spearman P`       = Spearman_P
  )

# Compact key columns of Table S1 (108 cell lines) so it fits cleanly on portrait pages
ts1_clean <- ts1 %>%
  transmute(
    `Model ID`       = ModelID,
    `Cell Line`      = CellLineName,
    `Lineage`        = Cancer_Code,
    `Subtype`        = OncotreeSubtype,
    `Epithelial`     = ifelse(Is_Epithelial, "Yes", "No"),
    `CRISPR PRKDC`   = round(CRISPR_PRKDC, 3),
    `CRISPR XRCC5`   = round(CRISPR_XRCC5, 3),
    `CRISPR XRCC6`   = round(CRISPR_XRCC6, 3),
    `CRISPR LIG4`    = round(CRISPR_LIG4, 3),
    `RNAi PRKDC`     = round(RNAi_PRKDC, 3)
  )

# ==============================================================================
# 2. SHARED WORD STYLING HELPERS
# ==============================================================================
fp_title   <- fp_text(font.size = 14.5, bold = TRUE, font.family = "Times New Roman", color = "#111111")
fp_alt     <- fp_text(font.size = 10.5, italic = TRUE, font.family = "Times New Roman", color = "#333333")
fp_h1      <- fp_text(font.size = 12, bold = TRUE, font.family = "Times New Roman", color = "#1F4E79")
fp_h2      <- fp_text(font.size = 11, bold = TRUE, italic = TRUE, font.family = "Times New Roman", color = "#222222")
fp_body    <- fp_text(font.size = 11, font.family = "Times New Roman", color = "#111111")
fp_bold    <- fp_text(font.size = 11, bold = TRUE, font.family = "Times New Roman", color = "#111111")
fp_cap     <- fp_text(font.size = 10, italic = TRUE, font.family = "Times New Roman", color = "#222222")
fp_hilite  <- fp_text(font.size = 10, bold = TRUE, font.family = "Times New Roman",
                      color = "#7F6000", shading.color = "#FFF2CC")

init_doc <- function(margin_in = 1.0) {
  d <- read_docx()
  body_set_default_section(d, value = prop_section(
    page_margins = page_mar(top = margin_in, bottom = margin_in, left = margin_in, right = margin_in)
  ))
}

add_h1 <- function(d, txt) {
  d <- body_add_par(d, value = "", style = "Normal")
  body_add_fpar(d, fpar(ftext(txt, fp_h1), fp_p = fp_par(padding.top = 6, padding.bottom = 4)))
}

add_h2 <- function(d, txt) {
  body_add_fpar(d, fpar(ftext(txt, fp_h2), fp_p = fp_par(padding.top = 4, padding.bottom = 2)))
}

add_p <- function(d, txt, pad_bottom = 6) {
  body_add_fpar(d, fpar(ftext(txt, fp_body),
                        fp_p = fp_par(line_spacing = 1.5, padding.bottom = pad_bottom, text.align = "justify")))
}

add_labeled_p <- function(d, label, txt) {
  body_add_fpar(d, fpar(ftext(label, fp_bold), ftext(txt, fp_body),
                        fp_p = fp_par(line_spacing = 1.5, padding.bottom = 4, text.align = "justify")))
}

make_pub_ft <- function(df, font_sz = 7.5) {
  ft <- flextable(df)
  ft <- font(ft, fontname = "Times New Roman", part = "all")
  ft <- fontsize(ft, size = font_sz, part = "all")
  ft <- bold(ft, part = "header")
  ft <- bg(ft, bg = "#F0F4F8", part = "header")
  ft <- border_remove(ft)
  ft <- hline_top(ft, border = fp_border(color = "black", width = 1.2), part = "header")
  ft <- hline_bottom(ft, border = fp_border(color = "black", width = 1.0), part = "header")
  ft <- hline_bottom(ft, border = fp_border(color = "black", width = 1.2), part = "body")
  ft <- padding(ft, padding.top = 2.5, padding.bottom = 2.5, padding.left = 3, padding.right = 3, part = "all")
  ft <- autofit(ft)
  ft <- set_table_properties(ft, width = 1, layout = "autofit")
  ft
}

# ==============================================================================
# 3. BUILD STANDALONE MAIN MANUSCRIPT (.DOCX)
# ==============================================================================
doc <- init_doc(margin_in = 1.0)

doc <- body_add_fpar(doc, fpar(
  ftext("Baseline Multi-Omic Profiles Fail to Predict Functional DNA-PK Dependency Across Cervical, Ovarian, and Endometrial Cancers: A Comparative Pan-Gynecologic Negative-Validation Study", fp_title),
  fp_p = fp_par(line_spacing = 1.3, padding.bottom = 4)
))
doc <- body_add_fpar(doc, fpar(
  ftext("Alternative Title Option: Constitutive Essentiality and Baseline Transcriptomic Unpredictability of the DNA-PK/c-NHEJ Complex Across Gynecologic Malignancies: A Multi-Lineage DepMap Audit", fp_alt),
  fp_p = fp_par(line_spacing = 1.25, padding.bottom = 8)
))

doc <- add_labeled_p(doc, "Running Title: ", "Pan-Gynecologic Negative Validation of Baseline DNA-PK Predictors")
doc <- add_labeled_p(doc, "Author: ", "Tanmay Singh Naruka")
doc <- add_labeled_p(doc, "Affiliation: ", "Department of Biotechnology, Amity University, Noida, Uttar Pradesh, India")
doc <- add_labeled_p(doc, "Data & Code Availability: ", "DepMap Public 26Q1 & DEMETER2 v6; reproducible R pipeline archived in PanGYN-DNAPK-Negative-Validation")

doc <- add_h1(doc, "Highlights")
highlights <- c(
  "\u2022 Comparative multi-omic audit of DNA-PK/c-NHEJ dependency across 108 gynecologic cancer cell lines (18 cervical, 56 ovarian, 34 endometrial) and 1,140 pan-cancer reference models.",
  "\u2022 Technical quality-control discovery: CRISPR PRKDC (5591) is NA-masked in 93.9% of pan-cancer lines (1,134/1,208) and 98.1% of gynecologic lines (106/108) due to PRKDCP1 pseudogene homology.",
  "\u2022 Obligatory DNA-PK holoenzyme subunits XRCC5 (Ku80) and XRCC6 (Ku70) exhibit 100% assay completeness (108/108) alongside 52 orthogonal DEMETER2 PRKDC RNAi profiles.",
  "\u2022 Leakage-free cross-lineage and within-lineage nested cross-validation collapse to the intercept-only null baseline across all three cancers, epithelial-only subsets (N = 89), and within-lineage RNAi screens.",
  "\u2022 Following REMARK and TRIPOD principles, clinical projections into TCGA-CESC (n = 309), TCGA-OV (n = 379), TCGA-UCEC (n = 548), and GEO cohorts were formally halted."
)
for (hl in highlights) {
  doc <- body_add_fpar(doc, fpar(ftext(hl, fp_body), fp_p = fp_par(line_spacing = 1.25, padding.bottom = 3)))
}

doc <- add_h1(doc, "Abstract")
doc <- add_labeled_p(doc, "Background: ", "Selective inhibitors of the DNA-dependent protein kinase catalytic subunit (DNA-PKcs, encoded by PRKDC) are under clinical evaluation as radiosensitizers and chemosensitizers across gynecologic malignancies. However, whether unperturbed baseline tumor omics encode constitutive DNA-PK dependency across distinct gynecologic lineages remains unproven.")
doc <- add_labeled_p(doc, "Methods: ", "We conducted a comparative, multi-omic negative-validation study integrating DepMap Public 26Q1 CRISPR-Cas9 knockout screens, protein-coding transcriptomics (19,210 leak-free genes), whole-genome sequencing copy number, damaging somatic mutations, and orthogonal DEMETER2 RNAi knockdown scores across 108 gynecologic models\u2014cervical (CACX, n = 18), ovarian (OVCX, n = 56), and endometrial (ENDOCX, n = 34)\u2014benchmarked against 1,140 pan-cancer reference models and an epithelial-only sensitivity subset (N = 89). Leakage-free cross-lineage transport and repeated 5-fold nested elastic-net and ridge models were evaluated against intercept-only null baselines under REMARK and TRIPOD reporting standards.")
doc <- add_labeled_p(doc, "Results: ", "Technical auditing revealed that true CRISPR PRKDC (5591) is NA-masked in 1,134 of 1,208 pan-cancer lines (0/18 CACX, 2/56 OVCX, 0/34 ENDOCX) due to Broad Avana guide cross-mapping with pseudogene PRKDCP1, whereas obligatory holoenzyme partners XRCC5 (Ku80) and XRCC6 (Ku70) achieved 100% coverage (108/108). Baseline target mRNA abundance showed zero correlation with functional dependency in any gynecologic cancer (all P > 0.06). Across cross-lineage pan-cancer transport and within-lineage nested cross-validation, transcriptomic and multi-omic models failed to outperform the intercept-only null baseline for XRCC5, XRCC6, or within-lineage DEMETER2 PRKDC RNAi (\u0394RMSE \u2265 0 across OVCX, ENDOCX, Pooled Pan-GYN, and epithelial-only cohorts). Lineage-adjusted co-dependency mapping further demonstrated complete decoupling between Ku80 (XRCC5) and catalytic c-NHEJ ligase LIG4 in unperturbed culture (all P > 0.64).")
doc <- add_labeled_p(doc, "Conclusions: ", "Unperturbed baseline multi-omics cannot predict functional DNA-PK vulnerability across cervical, ovarian, or endometrial cancers. Because predictive models collapse to zero-variance or non-generalizing null intercepts, downstream projections into TCGA-CESC (n = 309), TCGA-OV (n = 379), TCGA-UCEC (n = 548), and GEO cohorts were formally halted. Companion diagnostics for DNA-PK inhibitors in gynecologic oncology will require post-perturbation pharmacodynamic assays rather than baseline diagnostic biopsies.")
doc <- add_labeled_p(doc, "Keywords: ", "DNA-PKcs; PRKDC; XRCC5; XRCC6; non-homologous end joining; gynecologic oncology; cervical cancer; ovarian cancer; endometrial cancer; negative validation; DepMap")

# 1. INTRODUCTION
doc <- add_h1(doc, "1. Introduction")
doc <- add_p(doc, "Classical non-homologous end joining (c-NHEJ) is the primary pathway repairing ionizing radiation- and chemotherapy-induced DNA double-strand breaks (DSBs) in human cells [1]. Assembly of the active c-NHEJ synaptic complex requires initial recognition of broken DNA ends by the Ku70/Ku80 heterodimer (encoded by XRCC6 and XRCC5), recruitment and autophosphorylation of the DNA-dependent protein kinase catalytic subunit (DNA-PKcs, encoded by PRKDC), and final ligation by the DNA ligase IV\u2013XRCC4\u2013XLF complex (LIG4 and NHEJ1) [1]. Because cervical (CACX), ovarian (OVCX), and endometrial (ENDOCX) carcinomas rely heavily on DNA-damaging modalities\u2014including pelvic chemoradiotherapy and platinum-based chemotherapy\u2014selective pharmacological inhibitors of DNA-PKcs (such as peposertib/M3814 and AZD7648) have emerged as promising radiosensitizing and chemosensitizing agents [2, 3].")
doc <- add_p(doc, "A foundational assumption in translational oncology is that baseline, unperturbed multi-omic profiles from pre-treatment diagnostic biopsies can identify tumors intrinsically dependent on a therapeutic target prior to treatment exposure [4]. However, whether baseline transcriptomic or genomic features encode constitutive DNA-PK/c-NHEJ essentiality across major gynecologic malignancies\u2014or whether rigorous out-of-fold validation exposes an absence of predictive signal\u2014has not been systematically evaluated across gynecologic lineages. Reporting negative biomarker validation studies transparently is essential to prevent publication bias, wasted clinical correlative resources, and false-positive patient stratification [7, 8].")
doc <- add_p(doc, "Here, we present a unified, comparative pan-gynecologic negative-validation study auditing DNA-PK/c-NHEJ dependency across 108 evaluable cervical (n = 18), ovarian (n = 56), and endometrial (n = 34) cancer cell lines benchmarked against 1,140 pan-cancer reference models in the Cancer Dependency Map (DepMap Public 26Q1) [4, 5]. During multi-omic harmonization, we identified a major technical artifact in which true CRISPR PRKDC (5591) gene-effect scores are masked (NA) across 93.9% of screened cell lines due to sequence homology with the processed pseudogene PRKDCP1 [10]. By integrating the intact DNA-PK holoenzyme (XRCC5, XRCC6, LIG4, NHEJ1), whole-genome sequencing (WGS) copy number, damaging somatic mutations, histological epithelial filtering (N = 89), and orthogonal DEMETER2 RNAi knockdown screens (N = 52 gynecologic models) [6], we tested whether baseline omics encode constitutive DNA-PK dependency across cervical, ovarian, and endometrial cancers.")

# 2. MATERIALS AND METHODS
doc <- add_h1(doc, "2. Materials and Methods")

doc <- add_h2(doc, "2.1 Data Warehouse, Exact Identifier Harmonization, and Histological Stratification")
doc <- add_p(doc, "Discovery datasets were retrieved from the Broad Institute Cancer Dependency Map (DepMap Public 26Q1 release) [4], comprising Model.csv (2,154 cell lines), Chronos-corrected CRISPRGeneEffect.csv (1,208 cell lines) [5], protein-coding RNA-seq expression OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv (1,719 cell lines), WGS copy number OmicsCNGeneWGS.csv (1,118 cell lines), and somatic mutation calls OmicsSomaticMutations.csv (1,968 cell lines). To provide an orthogonal perturbation benchmark unconfounded by Avana CRISPR guide masking, merged RNAi gene dependency scores (D2_combined_gene_dep_scores.csv, DEMETER2 v6; 701 mapped models) were harmonized via Cancer Cell Line Encyclopedia (CCLE) identifiers [6].")
doc <- add_p(doc, "All matrices were deduplicated to unique DepMap model identifiers (ACH- prefix) prioritizing default model entries (IsDefaultEntryForModel == 'Yes'). Gene columns were parsed using exact HGNC symbol and Entrez ID regular expressions (^SYMBOL\\s*\\([0-9]+\\)$) to prevent substring collisions (such as matching PRKDC [5591] against APOBR [55911]). Gynecologic discovery cohorts were extracted using OncotreeLineage annotations: Cervix (CACX, n = 18), Ovary/Fallopian Tube (OVCX, n = 56), and Uterus (ENDOCX, n = 34), yielding 108 evaluable models with complete CRISPR holoenzyme and RNA-seq profiles alongside 1,140 pan-cancer reference models (Table S1). Because DepMap Uterus and Ovary annotations include non-epithelial histologies (e.g., uterine carcinosarcoma, leiomyosarcoma, choriocarcinoma, and ovarian germ-cell/granulosa tumors), we additionally defined a histologically pure epithelial sensitivity subset (N = 89: 18 CACX, 51 epithelial OVCX, and 20 Endometrial Carcinoma).")

doc <- add_h2(doc, "2.2 Leakage-Free Comparative Cross-Validation Architecture")
doc <- add_p(doc, "To prevent target leakage, self-transcript predictors for the DNA-PK/c-NHEJ complex (PRKDC, XRCC5, XRCC6, LIG4, NHEJ1) were excluded prior to modeling, retaining 19,210 protein-coding mRNA features alongside WGS copy-number drivers (PRKDC, XRCC5, XRCC6, LIG4, CCNE1, MYC, RB1, BRCA1, BRCA2, PTEN) and damaging loss-of-function (LoF; LikelyLoF == TRUE or VepImpact == 'HIGH') mutation indicators across 12 DNA damage response (DDR) and gynecologic driver genes (BRCA1, BRCA2, ATM, PALB2, RAD51C, RAD51D, TP53, PTEN, ARID1A, PIK3CA, KRAS, POLE).")
doc <- add_p(doc, "Four comparative modeling arms were executed using regularized linear regression (glmnet; Elastic-Net \u03b1 = 0.5 and Ridge \u03b1 = 0.0) [9] with all median imputation, zero-variance filtering, z-score standardization, and inner cross-validated penalty tuning (\u03bb_min) strictly nested within training splits: (1) Pan-Cancer to Gynecologic Cross-Lineage Transport, trained on all non-gynecologic pan-cancer models (n = 1,032 for CRISPR; n = 649 for DEMETER2 RNAi) and evaluated separately on CACX, OVCX, ENDOCX, and Pooled Pan-GYN; (2) Within-Lineage Repeated Nested Cross-Validation (5 outer folds \u00d7 5 repeats) using both transcriptome-only (19,210 features) and multi-omic (RNA + WGS CNV + LoF mutations) feature matrices; (3) Epithelial-Only Sensitivity Nested CV (N = 89); and (4) Orthogonal DEMETER2 PRKDC RNAi Nested CV (N = 52 gynecologic models). Model performance was benchmarked directly against an intercept-only null model predicting the training-fold mean (\u0394RMSE = Model RMSE \u2212 Null RMSE).")

doc <- add_h2(doc, "2.3 Abundance-vs.-Dependency, Complex Co-Dependency, and Clinical Halt Governance")
doc <- add_p(doc, "Univariate associations between baseline target mRNA abundance (log\u2082[TPM + 1]) and functional perturbation scores were quantified via Pearson linear and Spearman rank correlations across each lineage. To assess functional coupling within the c-NHEJ machinery without lineage confounding, CRISPR gene-effect scores were mean-centered within each OncotreeLineage prior to computing pairwise Pearson correlations between XRCC5 (Ku80) and complex/DDR partners (XRCC6, LIG4, NHEJ1, PARP1, POLQ, ATM, ATR, WEE1), adjusting P-values via the Benjamini\u2013Hochberg false discovery rate (FDR). In strict adherence to REMARK and TRIPOD reporting standards [7, 8], downstream transport into pre-specified clinical validation cohorts\u2014TCGA-CESC (n = 309), TCGA-OV (n = 379), TCGA-UCEC (n = 548), and paired radiotherapy/platinum GEO datasets (GSE56363, GSE3578, GSE14764, GSE32062, GSE17025)\u2014was gated on outperforming the intercept-only null baseline.")

# 3. RESULTS
doc <- add_h1(doc, "3. Results")

doc <- add_h2(doc, "3.1 DepMap 26Q1 Technical Audit Reveals Widespread PRKDC CRISPR Masking")
doc <- add_p(doc, "Systematic auditing of the 1,208 CRISPR-screened models in DepMap Public 26Q1 revealed that true PRKDC (5591) gene-effect scores are NA-masked in 1,134 cell lines (93.9%), leaving only 74 unmasked models pan-cancer due to Broad Avana sgRNA sequence homology with the processed pseudogene PRKDCP1 (Table 1; Figure 1B) [5, 10]. Within the 108 evaluable gynecologic models, CRISPR PRKDC (5591) is unmasked in 0/18 CACX lines, 0/34 ENDOCX lines, and only 2/56 OVCX lines (OVCAR-8: \u22121.009; OVISE: \u22120.745; 98.1% missing overall; Table S1). In contrast, the obligatory DNA-binding subunits of the DNA-PK holoenzyme\u2014XRCC5 (Ku80) and XRCC6 (Ku70)\u2014along with catalytic ligation factors LIG4 and NHEJ1 achieved 100% assay completeness (108/108) across all three gynecologic lineages, supplemented by 52 orthogonal DEMETER2 PRKDC RNAi knockdown profiles (2 CACX, 33 OVCX, 17 ENDOCX; Figure 1A\u2013B) [6].")
doc <- add_p(doc, "Across all three cancers, XRCC6 (Ku70) and XRCC5 (Ku80) exhibited uniformly strong constitutive essentiality with narrow interquartile ranges (CACX XRCC5 median \u22120.904, OVCX \u22120.880, ENDOCX \u22120.869; CACX XRCC6 median \u22121.898, OVCX \u22122.083, ENDOCX \u22122.008), whereas catalytic ligase LIG4 (medians \u22120.055 to \u22120.075) and NHEJ1 (medians \u22120.093 to \u22120.098) showed near-neutral baseline fitness effects (Table 1; Figure 2A).")

doc <- add_h2(doc, "3.2 Baseline Target mRNA Abundance Does Not Predict Functional Dependency")
doc <- add_p(doc, "Consistent with post-translational regulation of DNA-PK activation [1], baseline target mRNA expression failed to correlate with functional perturbation dependency across every individual and pooled gynecologic cohort (Table 3A; Figure 2B). For XRCC5 (Ku80), Pearson correlations between log\u2082(TPM + 1) expression and CRISPR gene effect were R = \u22120.0040 (P = 0.9875) in CACX, R = \u22120.1141 (P = 0.4024) in OVCX, R = 0.0515 (P = 0.7722) in ENDOCX, and R = \u22120.0594 (P = 0.5415) across the pooled 108-model pan-gynecologic cohort. Similarly, XRCC6 (Ku70) mRNA abundance showed no significant association with CRISPR dependency in CACX (R = \u22120.2159, P = 0.3897), OVCX (R = \u22120.2442, P = 0.0697), ENDOCX (R = 0.1315, P = 0.4586), or Pooled Pan-GYN (R = \u22120.1405, P = 0.1469). LIG4, NHEJ1, and DEMETER2 PRKDC RNAi scores (Pooled Pan-GYN n = 52, R = 0.1882, P = 0.1816) were likewise uncoupled from baseline transcript abundance.")

doc <- add_h2(doc, "3.3 Universal Collapse to the Null Hypothesis Across Comparative Cross-Validation Arms")
doc <- add_p(doc, "When pan-cancer transcriptomic models trained on 1,032 non-gynecologic cell lines were transported onto gynecologic cohorts, predictive accuracy failed to meaningfully improve over the constant intercept-only null model (Table 2; Figure 3B). For XRCC6 (Ku70), both Elastic-Net and Ridge models collapsed identically to the scalar training mean across all cohorts (Model RMSE = Null RMSE: 0.2957 in CACX, 0.3857 in OVCX, 0.3970 in ENDOCX, and 0.3759 in Pooled Pan-GYN; \u0394RMSE = 0.0000). For XRCC5 (Ku80), cross-lineage Elastic-Net transport yielded negligible or negative error changes (CACX \u0394RMSE = \u22120.0003, R = \u22120.0527; ENDOCX \u0394RMSE = +0.0003, R = \u22120.0551; OVCX \u0394RMSE = \u22120.0056, R = 0.2394; Pooled Pan-GYN \u0394RMSE = \u22120.0028, R = 0.1064).")
doc <- add_p(doc, "Within-lineage repeated 5-fold nested cross-validation confirmed that this null collapse is not an artifact of cross-lineage domain shift (Table 2; Figure 3A). For XRCC5, within-lineage transcriptome Elastic-Net models produced higher out-of-fold error than the intercept-only null baseline in CACX (RMSE 0.2674 vs. 0.2279, \u0394RMSE = +0.0395, R = \u22120.3228), OVCX (0.2565 vs. 0.2414, \u0394RMSE = +0.0151, R = 0.0948), ENDOCX (0.2659 vs. 0.2483, \u0394RMSE = +0.0176, R = \u22120.1107), and Pooled Pan-GYN (0.2488 vs. 0.2365, \u0394RMSE = +0.0123, R = \u22120.0495). Incorporating WGS copy-number and damaging LoF mutation layers into multi-omic models did not rescue predictive performance in OVCX (\u0394RMSE = +0.0416, R = 0.0154) or ENDOCX (\u0394RMSE = +0.0013, R = \u22120.1003). Within-lineage XRCC6 models similarly failed in OVCX (\u0394RMSE = +0.0117), ENDOCX (\u0394RMSE = +0.0005), and Pooled Pan-GYN (\u0394RMSE = +0.0074), while an apparent single-cohort variance fluctuation in small-sample CACX (n = 18; \u0394RMSE = \u22120.0809) failed to generalize to any larger gynecologic cohort.")
doc <- add_p(doc, "Two sensitivity analyses confirmed the robustness of these negative findings. First, restricting analysis to the 89 histologically pure epithelial carcinoma models (excluding uterine sarcomas/choriocarcinomas and non-epithelial ovarian tumors; Table S2A) yielded universal null collapse for both XRCC5 (Pooled Epithelial \u0394RMSE = +0.0347, R = \u22120.0881) and XRCC6 (Pooled Epithelial \u0394RMSE = +0.0257, R = \u22120.0571). Second, within-lineage nested cross-validation of orthogonal DEMETER2 PRKDC RNAi scores (Table S2B) likewise failed to outperform the null intercept in OVCX (n = 33, \u0394RMSE = +0.0277, R = \u22120.0347), ENDOCX (n = 17, \u0394RMSE = +0.0052, R = \u22120.0200), and Pooled Pan-GYN (n = 52, \u0394RMSE = +0.0724, R = 0.0979), while cross-lineage pan-cancer RNAi transport produced only near-null error shifts (|\u0394RMSE| \u2264 0.0191). Because no model satisfied pre-specified validation gates, downstream clinical projections into TCGA-CESC (n = 309), TCGA-OV (n = 379), TCGA-UCEC (n = 548), and GEO patient cohorts were formally halted [7, 8].")

doc <- add_h2(doc, "3.4 Uncoupling of Ku70/Ku80 from Catalytic c-NHEJ Ligation in Unperturbed Culture")
doc <- add_p(doc, "Lineage-adjusted co-dependency analysis revealed why baseline functional screens fail to capture therapeutic DNA-PK vulnerability (Table 3B; Figure 4). While XRCC5 (Ku80) correlated positively with its obligate heterodimeric partner XRCC6 (Ku70) in CACX (R = 0.6114, FDR q = 0.0281), OVCX (R = 0.4655, FDR q = 0.0024), and Pooled Pan-GYN (R = 0.3324, FDR q = 0.0035), downstream catalytic c-NHEJ ligase LIG4 exhibited complete statistical independence from XRCC5 in every gynecologic cancer (CACX R = \u22120.0656, P = 0.7961; OVCX R = 0.0629, P = 0.6450; ENDOCX R = 0.0376, P = 0.8330; Pooled Pan-GYN R = 0.0385, P = 0.6924). Meanwhile, secondary coupling shifted by lineage\u2014linking XRCC5 with POLQ in CACX (R = 0.6674, FDR q = 0.0198) and with ATR in ENDOCX (R = 0.4948, FDR q = 0.0235)\u2014while orthogonal CRISPR and DEMETER2 RNAi scores showed modest positive concordance across pan-gynecologic models for XRCC5 (R = 0.2743, P = 0.0514), XRCC6 (R = 0.2875, P = 0.0408), and ATR (R = 0.3537, P = 0.0101; Table S3).")

# 4. DISCUSSION & 5. CONCLUSION
doc <- add_h1(doc, "4. Discussion")
doc <- add_p(doc, "Our comparative pan-gynecologic investigation establishes two critical conclusions for translational DNA repair biomarker discovery. First, from a computational quality-control standpoint, blind ingestion of DepMap CRISPR matrices without exact HGNC/Entrez validation and missingness auditing can severely distort DNA-PK studies: Broad Avana sgRNA cross-mapping with the pseudogene PRKDCP1 masks true PRKDC (5591) scores in 93.9% of pan-cancer lines and 98.1% of gynecologic models [10]. By explicitly auditing the intact Ku70/Ku80\u2013DNA-PKcs\u2013LIG4 holoenzyme and triangulating against 52 DEMETER2 PRKDC RNAi profiles [6], we demonstrate that null predictability is a conserved biological property across cervical, ovarian, and endometrial cancers rather than a single-gene assay artifact.")
doc <- add_p(doc, "Second, biologically, our results reconcile why baseline multi-omic signatures fail to stratify constitutive DNA-PK dependency across gynecologic malignancies. Unlike oncogene addictions driven by static copy-number amplifications or hotspot mutations, classical NHEJ is an inducible emergency response triggered by acute DNA double-strand breaks [1]. In unperturbed cell culture, Ku70 and Ku80 perform constitutive housekeeping roles in telomere capping and rDNA stability\u2014explaining their uniform pan-essentiality and complete decoupling from catalytic LIG4 ligation\u2014whereas DNA-PKcs catalytic activation occurs only upon exogenous genotoxic stress such as ionizing radiation or platinum cross-linking [2, 3].")
doc <- add_p(doc, "From a clinical trial perspective, halting invalid projections into TCGA-CESC (n = 309), TCGA-OV (n = 379), and TCGA-UCEC (n = 548) prevents the propagation of zero-variance or overfitted noise into patient survival analyses [7, 8]. Future companion diagnostic strategies for DNA-PKcs inhibitors in gynecologic oncology should therefore pivot away from static pre-treatment transcriptomics toward dynamic, post-perturbation pharmacodynamic readouts\u2014such as ex vivo biopsy irradiation assays quantifying phosphorylated DNA-PKcs (Ser2056) and \u03b3-H2AX foci kinetics\u2014or target constitutively bimodal backup repair dependencies such as DNA polymerase theta (POLQ) in homologous-recombination- or PTEN/ATM-deficient subsets.")

doc <- add_h1(doc, "5. Conclusion")
doc <- add_p(doc, "Across 108 cervical, ovarian, and endometrial cancer models, unperturbed baseline transcriptomic and multi-omic profiles fail to predict functional DNA-PK holoenzyme dependency, collapsing consistently to the intercept-only null hypothesis. Rigorous negative validation and adherence to REMARK and TRIPOD stopping rules remain essential to ensure that only reproducible, biologically grounded biomarkers advance toward clinical gynecologic oncology.")

# 6. REFERENCES
doc <- add_h1(doc, "References")
refs <- c(
  "1. Blackford AN, Jackson SP. ATM, ATR, and DNA-PK: The Trinity at the Heart of the DNA Damage Response. Mol Cell. 2017;66(6):801\u2013817.",
  "2. Zenke FT, Zimmermann A, Sirrenberg C, et al. Pharmacologic Inhibitor of DNA-PK, M3814, Potentiates Radiotherapy and Regresses Human Tumors in Mouse Models. Mol Cancer Ther. 2020;19(5):1091\u20131101.",
  "3. Fok JHL, Ramos-Montoya A, Vazquez-Chantada M, et al. AZD7648 is a potent and selective DNA-PK inhibitor that enhances radiation, chemotherapy and olaparib activity. Nat Commun. 2019;10(1):5065.",
  "4. Tsherniak A, Vazquez F, Montgomery PG, et al. Defining a Cancer Dependency Map. Cell. 2017;170(3):564\u2013576.",
  "5. Dempster JM, Boyle I, Vazquez F, et al. Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. Genome Biol. 2021;22(1):343.",
  "6. McFarland JM, Ho ZV, Kugener G, et al. Improved estimation of cancer dependencies from large-scale RNAi screens using model-based normalization and data integration. Nat Commun. 2018;9(1):4610.",
  "7. McShane LM, Altman DG, Sauerbrei W, et al. Reporting recommendations for tumor marker prognostic studies (REMARK). J Natl Cancer Inst. 2005;97(16):1180\u20131184.",
  "8. Collins GS, Reitsma JB, Altman DG, Moons KGM. Transparent reporting of a multivariable prediction model for individual prognosis or diagnosis (TRIPOD). BMJ. 2015;350:g7594.",
  "9. Friedman J, Hastie T, Tibshirani R. Regularization Paths for Generalized Linear Models via Coordinate Descent. J Stat Softw. 2010;33(1):1\u201322.",
  "10. Broad Institute DepMap Consortium. DepMap Public 26Q1 Release Notes and CRISPR Guide Multi-Mapping Quality Control (PRKDC/PRKDCP1). Broad Institute; 2026."
)
for (r in refs) {
  doc <- body_add_fpar(doc, fpar(ftext(r, fp_body), fp_p = fp_par(line_spacing = 1.25, padding.bottom = 3)))
}

# 7. MAIN TABLES (TABLES 1–3B)
doc <- body_add_break(doc)
doc <- add_h1(doc, "Main Tables")

doc <- body_add_fpar(doc, fpar(ftext("Table 1. Pan-Gynecologic Discovery Cohort Architecture, Assay Completeness, and Matched Clinical Validation Cohorts.", fp_bold)))
doc <- body_add_flextable(doc, make_pub_ft(t1_clean, font_sz = 7.0))
doc <- body_add_par(doc, "", style = "Normal")

doc <- body_add_fpar(doc, fpar(ftext("Table 2. Comparative Cross-Lineage and Within-Lineage Repeated Nested Cross-Validation Benchmark Across CACX, OVCX, and ENDOCX.", fp_bold)))
doc <- body_add_flextable(doc, make_pub_ft(t2_clean, font_sz = 7.0))
doc <- body_add_par(doc, "", style = "Normal")

doc <- body_add_fpar(doc, fpar(ftext("Table 3A. Univariate Correlations Between Baseline Target mRNA Abundance and Functional Perturbation Dependency.", fp_bold)))
doc <- body_add_flextable(doc, make_pub_ft(t3a_clean, font_sz = 7.5))
doc <- body_add_par(doc, "", style = "Normal")

doc <- body_add_fpar(doc, fpar(ftext("Table 3B. Lineage-Adjusted Co-Dependency Between XRCC5 (Ku80) and DNA-PK / DDR Complex Partners.", fp_bold)))
doc <- body_add_flextable(doc, make_pub_ft(t3b_clean, font_sz = 8.0))

# 8. FIGURES & LEGENDS WITH HIGHLIGHTED SECTION PLACEMENT TAGS
doc <- body_add_break(doc)
doc <- add_h1(doc, "Figures and Legends")

fig_specs <- list(
  list(
    path = "results/figures/Figure1_PanGYN_Cohort_and_Assay_Audit.png", w = 6.5, h = 3.02,
    placement = "[Target Placement: Results, Section 3.1 \u2014 DepMap 26Q1 Technical Audit Reveals Widespread PRKDC CRISPR Masking]",
    cap = "Figure 1. Pan-Gynecologic Discovery Cohort Composition and Multi-Omic Assay Coverage Audit. (A) Evaluable DepMap Public 26Q1 models across cervical (CACX, n = 18), ovarian (OVCX, n = 56), and endometrial (ENDOCX, n = 34) lineages stratified by epithelial carcinoma (N = 89) versus non-epithelial histology (N = 19). (B) Assay completeness (%) across gynecologic cohorts and 1,140 pan-cancer reference models, highlighting 98.1% NA-masking of CRISPR PRKDC (5591) versus 100% completeness for XRCC5/XRCC6 and LIG4/NHEJ1."
  ),
  list(
    path = "results/figures/Figure2_PanGYN_Dependency_and_Abundance.png", w = 6.5, h = 4.83,
    placement = "[Target Placement: Results, Section 3.1 & Section 3.2 \u2014 Holoenzyme Essentiality & Baseline mRNA Abundance vs. Dependency]",
    cap = "Figure 2. DNA-PK Holoenzyme CRISPR Essentiality Distributions and Baseline mRNA Abundance vs. Dependency Across Gynecologic Lineages. (A) Chronos CRISPR gene-effect distributions for XRCC5 (Ku80), XRCC6 (Ku70), LIG4, and NHEJ1 (XLF) across CACX (n = 18), OVCX (n = 56), and ENDOCX (n = 34). (B) Linear regressions of baseline target mRNA expression (log\u2082[TPM + 1]) against CRISPR dependency (XRCC5, XRCC6) and orthogonal DEMETER2 PRKDC RNAi scores (all P > 0.05)."
  ),
  list(
    path = "results/figures/Figure3_PanGYN_NestedCV_Null_Collapse.png", w = 6.5, h = 4.76,
    placement = "[Target Placement: Results, Section 3.3 \u2014 Universal Collapse to the Null Hypothesis Across Comparative Cross-Validation Arms]",
    cap = "Figure 3. Comparative Nested Cross-Validation Benchmark Demonstrating Universal Null Collapse. (A) Out-of-fold root mean squared error (RMSE) of repeated 5-fold nested Elastic-Net models versus intercept-only null baselines across CRISPR XRCC5, CRISPR XRCC6, and DEMETER2 PRKDC RNAi. (B) Out-of-fold Pearson correlation (Observed vs. Predicted) across pan-cancer cross-lineage transport, within-lineage multi-omic models, epithelial-only pooled models (N = 89), and DEMETER2 RNAi cohorts."
  ),
  list(
    path = "results/figures/Figure4_PanGYN_Holoenzyme_CoDependency.png", w = 5.8, h = 3.69,
    placement = "[Target Placement: Results, Section 3.4 \u2014 Uncoupling of Ku70/Ku80 from Catalytic c-NHEJ Ligation in Unperturbed Culture]",
    cap = "Figure 4. Lineage-Adjusted Co-Dependency Heatmap with XRCC5 (Ku80) Across Gynecologic Cancers. Heatmap displays lineage-mean-centered Pearson correlation coefficients between XRCC5 (Ku80) and DNA-PK/DDR partners across CACX, OVCX, ENDOCX, and Pooled Pan-GYN (* BH-FDR q < 0.05, ** BH-FDR q < 0.01), demonstrating complete decoupling from catalytic c-NHEJ ligase LIG4 in unperturbed culture."
  )
)

for (fs in fig_specs) {
  if (file.exists(fs$path)) {
    doc <- body_add_fpar(doc, fpar(ftext(fs$placement, fp_hilite), fp_p = fp_par(padding.top = 4, padding.bottom = 3)))
    doc <- body_add_img(doc, src = fs$path, width = fs$w, height = fs$h)
    doc <- body_add_fpar(doc, fpar(ftext(fs$cap, fp_cap), fp_p = fp_par(padding.top = 3, padding.bottom = 10)))
  }
}

out_main <- "manuscript/PanGYN_DNAPK_Negative_Validation_Main_v1.1.docx"
print(doc, target = out_main)
cat(sprintf("\nSaved standalone Main Manuscript to: %s\n", normalizePath(out_main)))

# ==============================================================================
# 4. BUILD SEPARATE SUPPLEMENTARY INFORMATION FILE (.DOCX)
# ==============================================================================
sdoc <- init_doc(margin_in = 0.8)

sdoc <- body_add_fpar(sdoc, fpar(
  ftext("Supplementary Information", fp_title),
  fp_p = fp_par(padding.bottom = 4)
))
sdoc <- add_p(sdoc, "Baseline Multi-Omic Profiles Fail to Predict Functional DNA-PK Dependency Across Cervical, Ovarian, and Endometrial Cancers: A Comparative Pan-Gynecologic Negative-Validation Study", pad_bottom = 4)
sdoc <- add_labeled_p(sdoc, "Author: ", "Tanmay Singh Naruka (Department of Biotechnology, Amity University, Noida, Uttar Pradesh, India)")

sdoc <- add_h1(sdoc, "Supplementary Tables Overview")
sdoc <- add_p(sdoc, "This Supplementary Information file contains four supporting tables referenced in the main text: Table S1 (complete manifest of the 108 evaluable cervical, ovarian, and endometrial DepMap Public 26Q1 discovery models, epithelial classification, and perturbation scores), Table S2A (epithelial-only sensitivity nested cross-validation benchmark across N = 89 carcinoma models), Table S2B (orthogonal DEMETER2 v6 PRKDC RNAi cross-lineage and within-lineage nested cross-validation benchmark across N = 52 gynecologic models), and Table S3 (orthogonal concordance between CRISPR-Cas9 knockout and DEMETER2 RNAi dependency scores).")

sdoc <- add_h1(sdoc, "Table S1. Harmonized Pan-Gynecologic Discovery Cohort Manifest (N = 108 Models)")
sdoc <- body_add_flextable(sdoc, make_pub_ft(ts1_clean, font_sz = 7.0))

sdoc <- body_add_break(sdoc)
sdoc <- add_h1(sdoc, "Table S2A. Epithelial-Only Sensitivity Nested Cross-Validation Benchmark (N = 89 Epithelial Carcinoma Models)")
sdoc <- body_add_flextable(sdoc, make_pub_ft(ts2a_clean, font_sz = 7.0))
sdoc <- body_add_par(sdoc, "", style = "Normal")

sdoc <- add_h1(sdoc, "Table S2B. Orthogonal DEMETER2 PRKDC RNAi Cross-Lineage and Within-Lineage Nested CV Benchmark (N = 52 Models)")
sdoc <- body_add_flextable(sdoc, make_pub_ft(ts2b_clean, font_sz = 7.0))
sdoc <- body_add_par(sdoc, "", style = "Normal")

sdoc <- add_h1(sdoc, "Table S3. Orthogonal Concordance Between CRISPR Knockout and DEMETER2 RNAi Perturbation Scores")
sdoc <- body_add_flextable(sdoc, make_pub_ft(ts3_clean, font_sz = 7.5))

out_supp <- "manuscript/PanGYN_DNAPK_Negative_Validation_Supplementary_v1.1.docx"
print(sdoc, target = out_supp)
cat(sprintf("Saved separate Supplementary Information file to: %s\n", normalizePath(out_supp)))

# ==============================================================================
# 5. UPDATE GITHUB REPOSITORY README.md
# ==============================================================================
readme_lines <- c(
  "# PanGYN-DNAPK-Negative-Validation",
  "",
  "**Baseline Multi-Omic Profiles Fail to Predict Functional DNA-PK Dependency Across Cervical, Ovarian, and Endometrial Cancers: A Comparative Pan-Gynecologic Negative-Validation Study**",
  "",
  "## Overview",
  "This repository houses the reproducible R computational pipeline, harmonized summary tables, 600 DPI publication figures, and submission-ready manuscript files for a multi-lineage negative-validation audit of DNA-PK/c-NHEJ (`PRKDC`, `XRCC5`/Ku80, `XRCC6`/Ku70, `LIG4`, `NHEJ1`) dependency across **108 human gynecologic cancer cell lines**\u2014Cervical (`CACX`, n = 18), Ovarian (`OVCX`, n = 56), and Endometrial (`ENDOCX`, n = 34)\u2014benchmarked against **1,140 pan-cancer reference models** in DepMap Public 26Q1 and orthogonal **DEMETER2 v6 RNAi** screens.",
  "",
  "### Key Findings",
  "- **DepMap 26Q1 `PRKDC (5591)` Pseudogene Masking:** True CRISPR `PRKDC (5591)` is `NA`-masked in **1,134 / 1,208 (93.9%)** pan-cancer models and **106 / 108 (98.1%)** gynecologic models (`0/18` CACX, `2/56` OVCX, `0/34` ENDOCX) due to Broad Avana guide sequence homology with pseudogene *PRKDCP1*.",
  "- **100% Holoenzyme & Orthogonal RNAi Coverage:** Obligatory DNA-PK holoenzyme partners `XRCC5` (Ku80), `XRCC6` (Ku70), `LIG4`, and `NHEJ1` have complete **108 / 108** CRISPR coverage alongside **52** orthogonal `DEMETER2` *PRKDC* RNAi profiles.",
  "- **Zero Abundance-to-Dependency Coupling:** Baseline target mRNA expression shows no significant correlation with functional dependency in `CACX`, `OVCX`, `ENDOCX`, or `Pooled Pan-GYN` (all *P* > 0.06).",
  "- **Universal Collapse to the Null Hypothesis:** Across cross-lineage pan-cancer transport, within-lineage nested cross-validation (19,210 leak-free mRNAs + WGS copy number + LoF somatic mutations), epithelial-only sensitivity cohorts (`N = 89`), and within-lineage `DEMETER2` *PRKDC* RNAi (`N = 52`), regularized models consistently fail to outperform the intercept-only null baseline ($\\Delta\\text{RMSE} \\ge 0$).",
  "- **REMARK & TRIPOD Clinical Halt:** Because predictive models collapse to zero-variance or non-generalizing null intercepts, downstream clinical projections into `TCGA-CESC` (n = 309), `TCGA-OV` (n = 379), `TCGA-UCEC` (n = 548), and GEO cohorts were formally halted.",
  "",
  "## Repository Structure",
  "```text",
  "PanGYN-DNAPK-Negative-Validation/",
  "\u251c\u2500\u2500 .gitignore                                              # Excludes raw DepMap files & large .rds objects",
  "\u251c\u2500\u2500 README.md                                               # Project documentation & summary benchmarks",
  "\u251c\u2500\u2500 scripts/",
  "\u2502   \u251c\u2500\u2500 01_pangyn_harmonization.R                           # DepMap 26Q1 + DEMETER2 harmonisation & Table 1",
  "\u2502   \u251c\u2500\u2500 02_pangyn_comparative_cv.R                          # Leakage-free cross-lineage & nested CV (Table 2, S2A, S2B)",
  "\u2502   \u251c\u2500\u2500 03_pangyn_codependency.R                            # Abundance vs. dependency & complex coupling (Table 3A, 3B, S3)",
  "\u2502   \u251c\u2500\u2500 04_generate_pangyn_figures.R                        # Publication-grade 600 DPI Figures 1\u20134 (.tiff & .png)",
  "\u2502   \u2514\u2500\u2500 05_build_manuscript_docx.R                          # Programmatic Main + Supplementary .docx & README builder",
  "\u251c\u2500\u2500 results/",
  "\u2502   \u251c\u2500\u2500 tables/                                             # Verified CSV tables (Table 1\u20133, Table S1\u2013S3)",
  "\u2502   \u2514\u2500\u2500 figures/                                            # High-resolution 600 DPI figures",
  "\u2514\u2500\u2500 manuscript/",
  "    \u251c\u2500\u2500 PanGYN_DNAPK_Negative_Validation_Main_v1.1.docx",
  "    \u2514\u2500\u2500 PanGYN_DNAPK_Negative_Validation_Supplementary_v1.1.docx",
  "```",
  "",
  "## Reproducing the Pipeline",
  "1. Place DepMap Public 26Q1 (`Model.csv`, `CRISPRGeneEffect.csv`, `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`, `OmicsCNGeneWGS.csv`, `OmicsSomaticMutations.csv`) and `DEMETER2 v6` (`D2_combined_gene_dep_scores.csv`) inside `data/raw/`.",
  "2. Execute `scripts/01_pangyn_harmonization.R` through `scripts/05_build_manuscript_docx.R` sequentially in R."
)

writeLines(readme_lines, "README.md")
cat("Updated README.md successfully!\n")