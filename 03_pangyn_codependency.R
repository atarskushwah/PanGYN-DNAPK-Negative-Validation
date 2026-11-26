# scripts/03_pangyn_codependency.R
# Comparative Target Abundance vs. Dependency, DNA-PK Holoenzyme Coupling & DEMETER2 Concordance

library(tidyverse)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

cat("=== PHASE 3: ABUNDANCE VS. DEPENDENCY & DNA-PK COMPLEX COUPLING ===\n")
pangyn    <- readRDS("data/processed/01_pangyn_harmonized_108_models.rds")
pancancer <- readRDS("data/processed/01_pancancer_reference_models.rds")

safe_cor_stats <- function(x, y) {
  ok <- complete.cases(x, y)
  n_ok <- sum(ok)
  if (n_ok < 4 || sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    return(tibble(N = n_ok, Pearson_R = NA_real_, CI_Low = NA_real_, CI_High = NA_real_,
                  Pearson_P = NA_real_, Spearman_Rho = NA_real_, Spearman_P = NA_real_))
  }
  ct <- cor.test(x[ok], y[ok], method = "pearson")
  st <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  tibble(
    N            = n_ok,
    Pearson_R    = round(unname(ct$estimate), 4),
    CI_Low       = round(ct$conf.int[1], 4),
    CI_High      = round(ct$conf.int[2], 4),
    Pearson_P    = signif(ct$p.value, 4),
    Spearman_Rho = round(unname(st$estimate), 4),
    Spearman_P   = signif(st$p.value, 4)
  )
}

# ==============================================================================
# 1. TABLE 3A: TARGET mRNA ABUNDANCE VS. FUNCTIONAL DEPENDENCY ACROSS LINEAGES
# ==============================================================================
pairs_spec <- tribble(
  ~Target_Label,              ~RNA_Col,     ~Dep_Col,
  "XRCC5 (Ku80) CRISPR",      "RNA_XRCC5",  "CRISPR_XRCC5",
  "XRCC6 (Ku70) CRISPR",      "RNA_XRCC6",  "CRISPR_XRCC6",
  "LIG4 CRISPR",              "RNA_LIG4",   "CRISPR_LIG4",
  "NHEJ1 (XLF) CRISPR",       "RNA_NHEJ1",  "CRISPR_NHEJ1",
  "PRKDC (DNA-PKcs) DEMETER2","RNA_PRKDC",  "RNAi_PRKDC"
)

cohorts_list <- list(
  "CACX (n=18)"             = pangyn %>% filter(Cancer_Code == "CACX"),
  "OVCX (n=56)"             = pangyn %>% filter(Cancer_Code == "OVCX"),
  "ENDOCX (n=34)"           = pangyn %>% filter(Cancer_Code == "ENDOCX"),
  "Pooled Pan-GYN (N=108)"  = pangyn,
  "Pan-Cancer Ref (N=1140)" = pancancer
)

table3a <- map_dfr(names(cohorts_list), function(cname) {
  df <- cohorts_list[[cname]]
  pmap_dfr(pairs_spec, function(Target_Label, RNA_Col, Dep_Col) {
    res <- safe_cor_stats(df[[RNA_Col]], df[[Dep_Col]])
    bind_cols(tibble(Cohort = cname, Target = Target_Label), res)
  })
})

write_csv(table3a, "results/tables/Table3A_Abundance_vs_Dependency_PanGYN.csv")
cat("\n=== TABLE 3A: BASELINE mRNA ABUNDANCE VS. FUNCTIONAL DEPENDENCY ===\n")
print(as.data.frame(table3a), row.names = FALSE)

# ==============================================================================
# 2. TABLE 3B: LINEAGE-ADJUSTED DNA-PK HOLOENZYME & DDR CO-DEPENDENCY NETWORK
# ==============================================================================
# Mean-centre CRISPR scores within each OncotreeLineage to prevent lineage confounding
ddr_cols <- c("CRISPR_XRCC5", "CRISPR_XRCC6", "CRISPR_LIG4", "CRISPR_NHEJ1",
              "CRISPR_PARP1", "CRISPR_POLQ", "CRISPR_ATM", "CRISPR_ATR", "CRISPR_WEE1")

pancancer_centered <- pancancer %>%
  group_by(OncotreeLineage) %>%
  mutate(across(all_of(ddr_cols), ~ .x - mean(.x, na.rm = TRUE), .names = "CENT_{.col}")) %>%
  ungroup()

pangyn_centered <- pancancer_centered %>%
  filter(OncotreeLineage %in% c("Cervix", "Ovary/Fallopian Tube", "Uterus")) %>%
  mutate(Cancer_Code = case_when(
    OncotreeLineage == "Cervix"               ~ "CACX",
    OncotreeLineage == "Ovary/Fallopian Tube" ~ "OVCX",
    OncotreeLineage == "Uterus"               ~ "ENDOCX"
  ))

partner_genes <- c("CRISPR_XRCC6", "CRISPR_LIG4", "CRISPR_NHEJ1",
                   "CRISPR_PARP1", "CRISPR_POLQ", "CRISPR_ATM", "CRISPR_ATR", "CRISPR_WEE1")

codep_cohorts <- list(
  "CACX"           = pangyn_centered %>% filter(Cancer_Code == "CACX"),
  "OVCX"           = pangyn_centered %>% filter(Cancer_Code == "OVCX"),
  "ENDOCX"         = pangyn_centered %>% filter(Cancer_Code == "ENDOCX"),
  "Pooled Pan-GYN" = pangyn_centered
)

table3b <- map_dfr(names(codep_cohorts), function(cname) {
  df <- codep_cohorts[[cname]]
  map_dfr(partner_genes, function(pg) {
    res <- safe_cor_stats(df[["CENT_CRISPR_XRCC5"]], df[[paste0("CENT_", pg)]])
    tibble(
      Cohort       = cname,
      Anchor_Gene  = "XRCC5 (Ku80)",
      Partner_Gene = sub("^CRISPR_", "", pg),
      N            = res$N,
      Pearson_R    = res$Pearson_R,
      P_Value      = res$Pearson_P
    )
  }) %>%
    mutate(FDR_Q = signif(p.adjust(P_Value, method = "BH"), 4))
})

write_csv(table3b, "results/tables/Table3B_DNAPK_Holoenzyme_CoDependency.csv")
cat("\n=== TABLE 3B: LINEAGE-ADJUSTED CO-DEPENDENCY WITH XRCC5 (Ku80) ===\n")
print(as.data.frame(table3b), row.names = FALSE)

# ==============================================================================
# 3. TABLE S3: ORTHOGONAL DEMETER2 RNAi VS. CRISPR CONCORDANCE IN PAN-GYN
# ==============================================================================
table_s3 <- map_dfr(c("XRCC5", "XRCC6", "PARP1", "POLQ", "ATR", "WEE1"), function(g) {
  res <- safe_cor_stats(pangyn[[paste0("CRISPR_", g)]], pangyn[[paste0("RNAi_", g)]])
  bind_cols(tibble(Gene = g, Comparison = "CRISPR vs. DEMETER2 RNAi (Pan-GYN)"), res)
})

write_csv(table_s3, "results/tables/TableS3_DEMETER2_RNAi_vs_CRISPR_Concordance.csv")
cat("\n=== TABLE S3: ORTHOGONAL DEMETER2 RNAi VS. CRISPR CONCORDANCE ===\n")
print(as.data.frame(table_s3), row.names = FALSE)