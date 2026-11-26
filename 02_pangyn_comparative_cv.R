# scripts/02_pangyn_comparative_cv.R
# Comparative Cross-Lineage & Within-Lineage Nested CV Benchmark (CACX vs OVCX vs ENDOCX)

library(tidyverse)
library(glmnet)
library(rsample)

set.seed(2026)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

cat("=== PHASE 2: COMPARATIVE PAN-GYNECOLOGIC NESTED CV BENCHMARK ===\n")
pangyn    <- readRDS("data/processed/01_pangyn_harmonized_108_models.rds")
pancancer <- readRDS("data/processed/01_pancancer_reference_models.rds")

# Strictly exclude target complex mRNAs to prevent self-expression leakage
leakage_cols <- c("RNA_PRKDC", "RNA_XRCC5", "RNA_XRCC6", "RNA_LIG4", "RNA_NHEJ1")
rna_cols     <- setdiff(grep("^RNA_", names(pangyn), value = TRUE), leakage_cols)
cat(sprintf("Retained %d leak-free protein-coding transcriptomic predictors.\n", length(rna_cols)))

# Helper: Standardize training matrix and apply parameters to test matrix without leakage
prep_train_test <- function(X_tr_raw, X_te_raw) {
  for (j in seq_len(ncol(X_tr_raw))) {
    if (anyNA(X_tr_raw[, j]) || anyNA(X_te_raw[, j])) {
      med_val <- median(X_tr_raw[, j], na.rm = TRUE)
      if (is.na(med_val)) med_val <- 0
      X_tr_raw[is.na(X_tr_raw[, j]), j] <- med_val
      X_te_raw[is.na(X_te_raw[, j]), j] <- med_val
    }
  }
  tr_sds   <- apply(X_tr_raw, 2, sd)
  keep_col <- !is.na(tr_sds) & (tr_sds > 1e-8)
  X_tr_sub <- X_tr_raw[, keep_col, drop = FALSE]
  X_te_sub <- X_te_raw[, keep_col, drop = FALSE]
  
  X_tr <- scale(X_tr_sub)
  X_te <- scale(X_te_sub,
                center = attr(X_tr, "scaled:center"),
                scale  = attr(X_tr, "scaled:scale"))
  list(X_tr = X_tr, X_te = X_te)
}

# ==============================================================================
# 1. ARM A: PAN-CANCER CROSS-LINEAGE TRANSPORT ONTO CACX, OVCX & ENDOCX
# ==============================================================================
cat("\n[1/4] Running Pan-Cancer Cross-Lineage Transport (Non-GYN -> CACX / OVCX / ENDOCX)...\n")

non_gyn <- pancancer %>%
  filter(!(OncotreeLineage %in% c("Cervix", "Ovary/Fallopian Tube", "Uterus")))

run_cross_lineage_transport <- function(outcome_col, alpha_val = 0.5, algo_label = "ElasticNet (a=0.5)") {
  tr_df <- non_gyn %>% filter(!is.na(.data[[outcome_col]]))
  te_df <- pangyn  %>% filter(!is.na(.data[[outcome_col]]))
  
  mats   <- prep_train_test(as.matrix(tr_df[, rna_cols]), as.matrix(te_df[, rna_cols]))
  y_tr   <- tr_df[[outcome_col]]
  cv_fit <- cv.glmnet(mats$X_tr, y_tr, alpha = alpha_val, nfolds = 5, standardize = FALSE)
  
  preds  <- as.numeric(predict(cv_fit, newx = mats$X_te, s = "lambda.min"))
  coefs  <- as.numeric(coef(cv_fit, s = "lambda.min"))[-1]
  n_act  <- sum(coefs != 0)
  
  eval_tbl <- tibble(
    Cancer_Code = te_df$Cancer_Code,
    Observed    = te_df[[outcome_col]],
    Predicted   = preds,
    Null_Pred   = mean(y_tr)
  )
  
  by_cancer <- eval_tbl %>%
    group_by(Cancer_Code) %>%
    summarise(
      N_Eval          = n(),
      Model_RMSE      = round(sqrt(mean((Observed - Predicted)^2)), 4),
      Null_RMSE       = round(sqrt(mean((Observed - Null_Pred)^2)), 4),
      Delta_RMSE      = round(Model_RMSE - Null_RMSE, 4),
      Transport_R     = round(suppressWarnings(cor(Observed, Predicted)), 4),
      Active_Features = n_act,
      .groups         = "drop"
    )
  
  pooled_row <- eval_tbl %>%
    summarise(
      Cancer_Code     = "Pooled Pan-GYN",
      N_Eval          = n(),
      Model_RMSE      = round(sqrt(mean((Observed - Predicted)^2)), 4),
      Null_RMSE       = round(sqrt(mean((Observed - Null_Pred)^2)), 4),
      Delta_RMSE      = round(Model_RMSE - Null_RMSE, 4),
      Transport_R     = round(suppressWarnings(cor(Observed, Predicted)), 4),
      Active_Features = n_act
    )
  
  bind_rows(by_cancer, pooled_row) %>%
    mutate(
      Validation_Arm = "Pan-Cancer -> GYN Transport",
      Outcome        = outcome_col,
      Algorithm      = algo_label,
      .before        = 1
    )
}

transport_res <- bind_rows(
  run_cross_lineage_transport("CRISPR_XRCC5", 0.5, "Transcriptome ElasticNet (a=0.5)"),
  run_cross_lineage_transport("CRISPR_XRCC5", 0.0, "Transcriptome Ridge (a=0.0)"),
  run_cross_lineage_transport("CRISPR_XRCC6", 0.5, "Transcriptome ElasticNet (a=0.5)"),
  run_cross_lineage_transport("CRISPR_XRCC6", 0.0, "Transcriptome Ridge (a=0.0)"),
  run_cross_lineage_transport("RNAi_PRKDC",   0.5, "Transcriptome ElasticNet (a=0.5)"),
  run_cross_lineage_transport("RNAi_PRKDC",   0.0, "Transcriptome Ridge (a=0.0)")
)

# ==============================================================================
# 2. ARM B: WITHIN-LINEAGE REPEATED 5-FOLD NESTED CV (5 FOLDS x 5 REPEATS)
# ==============================================================================
cat("[2/4] Running Within-Lineage Repeated Nested CV across CACX, OVCX, ENDOCX & Pooled Pan-GYN...\n")

run_within_nested_cv <- function(df, outcome_col, feat_cols, cohort_label,
                                 algo_label = "Transcriptome ElasticNet (a=0.5)", alpha_val = 0.5) {
  sub_df <- df %>% filter(!is.na(.data[[outcome_col]]))
  X_all  <- as.matrix(sub_df[, feat_cols, drop = FALSE])
  y_vec  <- sub_df[[outcome_col]]
  
  n_inner <- min(5, max(3, floor(nrow(sub_df) * 0.8 / 3)))
  folds   <- vfold_cv(sub_df, v = 5, repeats = 5)
  
  oof_tbl <- map_dfr(seq_len(nrow(folds)), function(i) {
    tr_idx <- as.integer(folds$splits[[i]])
    te_idx <- setdiff(seq_len(nrow(sub_df)), tr_idx)
    
    mats   <- prep_train_test(X_all[tr_idx, , drop = FALSE], X_all[te_idx, , drop = FALSE])
    y_tr   <- y_vec[tr_idx]
    y_te   <- y_vec[te_idx]
    
    cv_fit <- cv.glmnet(mats$X_tr, y_tr, alpha = alpha_val, nfolds = n_inner, standardize = FALSE)
    preds  <- as.numeric(predict(cv_fit, newx = mats$X_te, s = "lambda.min"))
    coefs  <- as.numeric(coef(cv_fit, s = "lambda.min"))[-1]
    
    tibble(
      Repeat    = folds$id[i],
      Observed  = y_te,
      Predicted = preds,
      Null_Pred = rep(mean(y_tr), length(y_te)),
      N_Active  = sum(coefs != 0)
    )
  })
  
  rep_summary <- oof_tbl %>%
    group_by(Repeat) %>%
    summarise(
      RMSE_Model = sqrt(mean((Observed - Predicted)^2)),
      RMSE_Null  = sqrt(mean((Observed - Null_Pred)^2)),
      OOF_R      = suppressWarnings(cor(Observed, Predicted)),
      Med_Act    = median(N_Active),
      .groups    = "drop"
    )
  
  tibble(
    Validation_Arm  = "Within-Lineage Nested CV (5x5)",
    Outcome         = outcome_col,
    Algorithm       = algo_label,
    Cancer_Code     = cohort_label,
    N_Eval          = nrow(sub_df),
    Model_RMSE      = round(mean(rep_summary$RMSE_Model), 4),
    Null_RMSE       = round(mean(rep_summary$RMSE_Null), 4),
    Delta_RMSE      = round(Model_RMSE - Null_RMSE, 4),
    Transport_R     = round(mean(rep_summary$OOF_R, na.rm = TRUE), 4),
    Active_Features = median(rep_summary$Med_Act)
  )
}

multiomic_cols <- c(rna_cols, grep("^(CN_|MUT_)", names(pangyn), value = TRUE))

within_res <- bind_rows(
  # XRCC5 (Ku80) across CACX, OVCX, ENDOCX, Pooled Pan-GYN
  run_within_nested_cv(pangyn %>% filter(Cancer_Code == "CACX"),   "CRISPR_XRCC5", rna_cols, "CACX",           "Transcriptome ElasticNet (a=0.5)", 0.5),
  run_within_nested_cv(pangyn %>% filter(Cancer_Code == "OVCX"),   "CRISPR_XRCC5", rna_cols, "OVCX",           "Transcriptome ElasticNet (a=0.5)", 0.5),
  run_within_nested_cv(pangyn %>% filter(Cancer_Code == "ENDOCX"), "CRISPR_XRCC5", rna_cols, "ENDOCX",         "Transcriptome ElasticNet (a=0.5)", 0.5),
  run_within_nested_cv(pangyn,                                     "CRISPR_XRCC5", rna_cols, "Pooled Pan-GYN", "Transcriptome ElasticNet (a=0.5)", 0.5),
  run_within_nested_cv(pangyn %>% filter(Cancer_Code == "OVCX"),   "CRISPR_XRCC5", multiomic_cols, "OVCX",     "Multi-Omic (RNA+CNV+MUT) ENet",    0.5),
  run_within_nested_cv(pangyn %>% filter(Cancer_Code == "ENDOCX"), "CRISPR_XRCC5", multiomic_cols, "ENDOCX",   "Multi-Omic (RNA+CNV+MUT) ENet",    0.5),
  
  # XRCC6 (Ku70) across CACX, OVCX, ENDOCX, Pooled Pan-GYN
  run_within_nested_cv(pangyn %>% filter(Cancer_Code == "CACX"),   "CRISPR_XRCC6", rna_cols, "CACX",           "Transcriptome ElasticNet (a=0.5)", 0.5),
  run_within_nested_cv(pangyn %>% filter(Cancer_Code == "OVCX"),   "CRISPR_XRCC6", rna_cols, "OVCX",           "Transcriptome ElasticNet (a=0.5)", 0.5),
  run_within_nested_cv(pangyn %>% filter(Cancer_Code == "ENDOCX"), "CRISPR_XRCC6", rna_cols, "ENDOCX",         "Transcriptome ElasticNet (a=0.5)", 0.5),
  run_within_nested_cv(pangyn,                                     "CRISPR_XRCC6", rna_cols, "Pooled Pan-GYN", "Transcriptome ElasticNet (a=0.5)", 0.5)
)

table2 <- bind_rows(transport_res %>% filter(Outcome != "RNAi_PRKDC"), within_res)
write_csv(table2, "results/tables/Table2_PanGYN_Nested_CV_Benchmark.csv")

# ==============================================================================
# 3. ARM C: EPITHELIAL-ONLY SENSITIVITY BENCHMARK (N = 89) -> TABLE S2A
# ==============================================================================
cat("[3/4] Running Epithelial-Only Sensitivity Benchmark (N = 89)...\n")
epi_df <- pangyn %>% filter(Is_Epithelial)

table_s2a <- bind_rows(
  run_within_nested_cv(epi_df %>% filter(Cancer_Code == "OVCX"),   "CRISPR_XRCC5", rna_cols, "Epithelial OVCX (n=51)",     "Transcriptome ElasticNet (a=0.5)", 0.5),
  run_within_nested_cv(epi_df %>% filter(Cancer_Code == "ENDOCX"), "CRISPR_XRCC5", rna_cols, "Endometrial Carc. (n=20)",   "Transcriptome ElasticNet (a=0.5)", 0.5),
  run_within_nested_cv(epi_df,                                     "CRISPR_XRCC5", rna_cols, "Pooled Epithelial GYN (89)", "Transcriptome ElasticNet (a=0.5)", 0.5),
  run_within_nested_cv(epi_df %>% filter(Cancer_Code == "OVCX"),   "CRISPR_XRCC6", rna_cols, "Epithelial OVCX (n=51)",     "Transcriptome ElasticNet (a=0.5)", 0.5),
  run_within_nested_cv(epi_df %>% filter(Cancer_Code == "ENDOCX"), "CRISPR_XRCC6", rna_cols, "Endometrial Carc. (n=20)",   "Transcriptome ElasticNet (a=0.5)", 0.5),
  run_within_nested_cv(epi_df,                                     "CRISPR_XRCC6", rna_cols, "Pooled Epithelial GYN (89)", "Transcriptome ElasticNet (a=0.5)", 0.5)
)
write_csv(table_s2a, "results/tables/TableS2A_Epithelial_Sensitivity_CV.csv")

# ==============================================================================
# 4. ARM D: ORTHOGONAL DEMETER2 PRKDC RNAi BENCHMARK (N = 52 GYN) -> TABLE S2B
# ==============================================================================
cat("[4/4] Running Orthogonal DEMETER2 PRKDC RNAi Benchmark (OVCX n=33, ENDOCX n=17, Pooled n=52)...\n")

table_s2b <- bind_rows(
  transport_res %>% filter(Outcome == "RNAi_PRKDC"),
  run_within_nested_cv(pangyn %>% filter(Cancer_Code == "OVCX"),   "RNAi_PRKDC", rna_cols, "OVCX RNAi (n=33)",       "Transcriptome ElasticNet (a=0.5)", 0.5),
  run_within_nested_cv(pangyn %>% filter(Cancer_Code == "ENDOCX"), "RNAi_PRKDC", rna_cols, "ENDOCX RNAi (n=17)",     "Transcriptome ElasticNet (a=0.5)", 0.5),
  run_within_nested_cv(pangyn,                                     "RNAi_PRKDC", rna_cols, "Pooled GYN RNAi (n=52)", "Transcriptome ElasticNet (a=0.5)", 0.5)
)
write_csv(table_s2b, "results/tables/TableS2B_DEMETER2_PRKDC_RNAi_CV.csv")

cat("\n=== TABLE 2: COMPARATIVE CROSS-LINEAGE & WITHIN-LINEAGE NESTED CV ===\n")
print(as.data.frame(table2), row.names = FALSE)

cat("\n=== TABLE S2A: EPITHELIAL-ONLY SENSITIVITY BENCHMARK ===\n")
print(as.data.frame(table_s2a), row.names = FALSE)

cat("\n=== TABLE S2B: ORTHOGONAL DEMETER2 PRKDC RNAi BENCHMARK ===\n")
print(as.data.frame(table_s2b), row.names = FALSE)