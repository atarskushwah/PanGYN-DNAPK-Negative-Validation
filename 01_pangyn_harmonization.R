# scripts/01_pangyn_harmonization.R
# Unified Pan-Gynecologic (CACX vs. OVCX vs. ENDOCX) DepMap 26Q1 + DEMETER2 Harmonisation

library(tidyverse)
library(data.table)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# 1. Create .gitignore immediately so GitHub pushes stay well under 100 MB
writeLines(c(
  "# R & macOS artifacts",
  ".Rhistory",
  ".RData",
  ".Rproj.user/",
  ".DS_Store",
  "",
  "# Large raw and processed omics files",
  "data/raw/",
  "data/processed/*.rds",
  "*.zip",
  "*.gz"
), ".gitignore")

# ==============================================================================
# 2. EXACT & DEDUPLICATED DEPMAP PUBLIC 26Q1 LOADER
# ==============================================================================
load_depmap_exact <- function(path, label, deduplicate = TRUE) {
  df <- fread(path, data.table = FALSE, check.names = FALSE)
  
  if ("DepMap_ID" %in% names(df) && !("ModelID" %in% names(df))) {
    names(df)[names(df) == "DepMap_ID"] <- "ModelID"
  } else if (!("ModelID" %in% names(df))) {
    names(df)[1] <- "ModelID"
  }
  
  id_idx <- which(names(df) == "ModelID")
  if (length(id_idx) > 1) df <- df[, -id_idx[-1], drop = FALSE]
  
  df <- df[!is.na(df$ModelID) & grepl("^ACH-", df$ModelID), , drop = FALSE]
  
  if (deduplicate) {
    if ("IsDefaultEntryForModel" %in% names(df)) {
      is_def <- as.character(df$IsDefaultEntryForModel) %in% c("Yes", "TRUE", "True", "1")
      df <- df[order(!is_def), , drop = FALSE]
    }
    df <- df[!duplicated(df$ModelID), , drop = FALSE]
  }
  
  cat(sprintf("[%-10s] %d unique ACH- models loaded\n", label, length(unique(df$ModelID))))
  return(df)
}

cat("=== PHASE 1: PAN-GYNECOLOGIC MULTI-OMIC HARMONISATION ===\n")
models <- load_depmap_exact("data/raw/Model.csv", "Model")
crispr <- load_depmap_exact("data/raw/CRISPRGeneEffect.csv", "CRISPR")
expr   <- load_depmap_exact("data/raw/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", "Expression")
cnv    <- load_depmap_exact("data/raw/OmicsCNGeneWGS.csv", "CNV")
muts   <- load_depmap_exact("data/raw/OmicsSomaticMutations.csv", "Mutations", deduplicate = FALSE)

get_exact_col <- function(df, symbol) {
  hit <- grep(paste0("^", symbol, "\\s*\\([0-9]+\\)$"), names(df), value = TRUE)
  if (length(hit) == 0) return(NA_character_)
  return(hit[1])
}

# ==============================================================================
# 3. LOAD & TRANSPOSE DEMETER2 COMBINED RNAi (D2_combined_gene_dep_scores.csv)
# ==============================================================================
dnapk_genes <- c("PRKDC", "XRCC5", "XRCC6", "LIG4", "NHEJ1", "PARP1", "POLQ", "ATM", "ATR", "WEE1")
rnai_path   <- "data/raw/D2_combined_gene_dep_scores.csv"

if (file.exists(rnai_path)) {
  d2_raw <- fread(rnai_path, data.table = FALSE, check.names = FALSE)
  names(d2_raw)[1] <- "GeneRow"
  
  # Extract target DDR gene rows (handles both "PRKDC (5591)" and readthrough symbols)
  d2_hits <- map_dfr(dnapk_genes, function(g) {
    idx <- grep(paste0("(^|&)", g, "\\s*\\("), d2_raw$GeneRow)
    if (length(idx) == 0) return(NULL)
    row_vals <- as.numeric(d2_raw[idx[1], -1])
    tibble(
      CCLEName = names(d2_raw)[-1],
      Gene     = paste0("RNAi_", g),
      Score    = row_vals
    )
  })
  
  # Map CCLEName columns in DEMETER2 to DepMap ModelID via Model.csv
  ccle_lookup <- models %>%
    filter(!is.na(CCLEName) & CCLEName != "") %>%
    distinct(CCLEName, .keep_all = TRUE) %>%
    select(CCLEName, ModelID)
  
  rnai_sub <- d2_hits %>%
    inner_join(ccle_lookup, by = "CCLEName") %>%
    select(ModelID, Gene, Score) %>%
    pivot_wider(names_from = Gene, values_from = Score)
  
  cat(sprintf("[%-10s] %d unique ACH- models mapped from DEMETER2 RNAi\n", "DEMETER2", nrow(rnai_sub)))
} else {
  warning("D2_combined_gene_dep_scores.csv not found in data/raw/; proceeding without RNAi layer.")
  rnai_sub <- tibble(ModelID = character())
}

# ==============================================================================
# 4. DEFINE GYNECOLOGIC COHORTS & EPITHELIAL-ONLY SENSITIVITY FLAGS
# ==============================================================================
gyn_meta <- models %>%
  filter(OncotreeLineage %in% c("Cervix", "Ovary/Fallopian Tube", "Uterus")) %>%
  mutate(
    Cancer_Code = case_when(
      OncotreeLineage == "Cervix"               ~ "CACX",
      OncotreeLineage == "Ovary/Fallopian Tube" ~ "OVCX",
      OncotreeLineage == "Uterus"               ~ "ENDOCX"
    ),
    Is_Epithelial = case_when(
      Cancer_Code == "CACX" ~ TRUE,
      Cancer_Code == "OVCX" & !grepl("Germ Cell|Small Cell|Granulosa|Brenner|Sarcoma", OncotreeSubtype, ignore.case = TRUE) ~ TRUE,
      Cancer_Code == "ENDOCX" & OncotreeSubtype == "Endometrial Carcinoma" ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  select(ModelID, CellLineName, StrippedCellLineName, CCLEName, Cancer_Code,
         OncotreeLineage, OncotreePrimaryDisease, OncotreeSubtype, Is_Epithelial)

# ==============================================================================
# 5. EXTRACT CRISPR, TRANSCRIPTOME, WGS COPY NUMBER & LOF MUTATIONS
# ==============================================================================
crispr_map <- setNames(sapply(dnapk_genes, function(g) get_exact_col(crispr, g)), dnapk_genes)
crispr_map <- crispr_map[!is.na(crispr_map)]
crispr_sub <- crispr[, c("ModelID", unname(crispr_map)), drop = FALSE]
names(crispr_sub)[-1] <- paste0("CRISPR_", names(crispr_map))
crispr_sub[-1] <- lapply(crispr_sub[-1], as.numeric)

# Clean protein-coding transcriptome
rna_gene_cols <- grep("\\([0-9]+\\)$", names(expr), value = TRUE)
expr_mat <- as.matrix(expr[, rna_gene_cols])
rownames(expr_mat) <- expr$ModelID
clean_rna_sym <- gsub("\\s*\\([0-9]+\\)$", "", colnames(expr_mat))
keep_rna <- !duplicated(clean_rna_sym)
expr_mat <- expr_mat[, keep_rna, drop = FALSE]
colnames(expr_mat) <- paste0("RNA_", clean_rna_sym[keep_rna])
expr_df <- as_tibble(expr_mat, rownames = "ModelID")

# WGS Copy Number for DNA-PK & GYN drivers
cn_genes <- c("PRKDC", "XRCC5", "XRCC6", "LIG4", "CCNE1", "MYC", "RB1", "BRCA1", "BRCA2", "PTEN")
cn_map   <- setNames(sapply(cn_genes, function(g) get_exact_col(cnv, g)), cn_genes)
cn_map   <- cn_map[!is.na(cn_map)]
cnv_sub  <- cnv[, c("ModelID", unname(cn_map)), drop = FALSE]
names(cnv_sub)[-1] <- paste0("CN_", names(cn_map))
cnv_sub[-1] <- lapply(cnv_sub[-1], as.numeric)

# Damaging Loss-of-Function somatic mutations
gene_col_mut  <- intersect(c("HugoSymbol", "Hugo_Symbol"), names(muts))[1]
ddr_mut_genes <- c("BRCA1", "BRCA2", "ATM", "PALB2", "RAD51C", "RAD51D",
                   "TP53", "PTEN", "ARID1A", "PIK3CA", "KRAS", "POLE")

mut_sub <- muts[muts$ModelID %in% gyn_meta$ModelID & muts[[gene_col_mut]] %in% ddr_mut_genes, , drop = FALSE]
if ("LikelyLoF" %in% names(mut_sub)) {
  mut_sub <- mut_sub[as.character(mut_sub$LikelyLoF) %in% c("TRUE", "True", "1"), , drop = FALSE]
} else if ("VepImpact" %in% names(mut_sub)) {
  mut_sub <- mut_sub[mut_sub$VepImpact == "HIGH", , drop = FALSE]
}

mut_matrix <- mut_sub %>%
  transmute(ModelID, Gene = .data[[gene_col_mut]], Status = 1L) %>%
  distinct() %>%
  pivot_wider(names_from = Gene, values_from = Status, values_fill = 0L, names_prefix = "MUT_")

# ==============================================================================
# 6. MERGE DISCOVERY COHORTS (PAN-GYN N = 108 & PAN-CANCER REFERENCE)
# ==============================================================================
pangyn_cohort <- gyn_meta %>%
  inner_join(crispr_sub, by = "ModelID") %>%
  inner_join(expr_df,    by = "ModelID") %>%
  left_join(rnai_sub,    by = "ModelID") %>%
  left_join(cnv_sub,     by = "ModelID") %>%
  left_join(mut_matrix,  by = "ModelID") %>%
  mutate(across(starts_with("MUT_"), ~replace_na(., 0L))) %>%
  filter(!is.na(CRISPR_XRCC5) & !is.na(CRISPR_XRCC6))

pancancer_cohort <- models %>%
  select(ModelID, CellLineName, OncotreeLineage, OncotreeSubtype) %>%
  inner_join(crispr_sub, by = "ModelID") %>%
  inner_join(expr_df,    by = "ModelID") %>%
  left_join(rnai_sub,    by = "ModelID") %>%
  filter(!is.na(CRISPR_XRCC5) & !is.na(CRISPR_XRCC6))

saveRDS(pangyn_cohort,    "data/processed/01_pangyn_harmonized_108_models.rds")
saveRDS(pancancer_cohort, "data/processed/01_pancancer_reference_models.rds")

# ==============================================================================
# 7. GENERATE TABLE 1 & CELL-LINE MANIFEST (TABLE S1)
# ==============================================================================
table1 <- pangyn_cohort %>%
  group_by(Cancer_Code, OncotreeLineage) %>%
  summarise(
    Evaluable_CRISPR_RNA    = n(),
    Epithelial_Only_N       = sum(Is_Epithelial),
    CRISPR_PRKDC_Unmasked_N = sum(!is.na(CRISPR_PRKDC)),
    DEMETER2_PRKDC_RNAi_N   = sum(!is.na(RNAi_PRKDC)),
    XRCC5_Ku80_Median_IQR   = sprintf("%.3f (%.3f to %.3f)",
                                      median(CRISPR_XRCC5),
                                      quantile(CRISPR_XRCC5, 0.25),
                                      quantile(CRISPR_XRCC5, 0.75)),
    XRCC6_Ku70_Median_IQR   = sprintf("%.3f (%.3f to %.3f)",
                                      median(CRISPR_XRCC6),
                                      quantile(CRISPR_XRCC6, 0.25),
                                      quantile(CRISPR_XRCC6, 0.75)),
    WGS_CNV_Overlap_N       = sum(!is.na(CN_XRCC5)),
    .groups                 = "drop"
  ) %>%
  mutate(
    Matched_TCGA_Cohort = case_when(
      Cancer_Code == "CACX"   ~ "TCGA-CESC (n = 309 primary tumours)",
      Cancer_Code == "OVCX"   ~ "TCGA-OV (n = 379 primary HGSOC tumours)",
      Cancer_Code == "ENDOCX" ~ "TCGA-UCEC (n = 548 primary uterine tumours)"
    ),
    Matched_GEO_Cohorts = case_when(
      Cancer_Code == "CACX"   ~ "GSE56363 (n = 21 CRT), GSE3578 (39 paired biopsies)",
      Cancer_Code == "OVCX"   ~ "GSE14764 / GSE32062 (platinum-treated HGSOC)",
      Cancer_Code == "ENDOCX" ~ "GSE17025 / TCGA-UCEC molecular subtypes (POLE/MSI/CN)"
    )
  )

write_csv(table1, "results/tables/Table1_PanGYN_Cohort_Architecture.csv")
write_csv(pangyn_cohort %>% select(ModelID, CellLineName, Cancer_Code, OncotreeSubtype, Is_Epithelial,
                                   starts_with("CRISPR_"), starts_with("RNAi_"),
                                   starts_with("CN_"), starts_with("MUT_")),
          "results/tables/TableS1_PanGYN_108_CellLines_Manifest.csv")

cat("\n=== TABLE 1: PAN-GYNECOLOGIC DISCOVERY & VALIDATION ARCHITECTURE ===\n")
print(as.data.frame(table1), row.names = FALSE)
cat(sprintf("\nPan-Cancer Reference Models: %d | Pan-GYN Models: %d (Epithelial-Only: %d)\n",
            nrow(pancancer_cohort), nrow(pangyn_cohort), sum(pangyn_cohort$Is_Epithelial)))