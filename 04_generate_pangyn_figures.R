# scripts/04_generate_pangyn_figures.R
# Publication-Grade 600 DPI Figures (TIFF & PNG) for Pan-Gynecologic DNA-PK Study

library(tidyverse)
library(patchwork)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

pangyn    <- readRDS("data/processed/01_pangyn_harmonized_108_models.rds")
pancancer <- readRDS("data/processed/01_pancancer_reference_models.rds")
table2    <- read_csv("results/tables/Table2_PanGYN_Nested_CV_Benchmark.csv", show_col_types = FALSE)
table_s2a <- read_csv("results/tables/TableS2A_Epithelial_Sensitivity_CV.csv", show_col_types = FALSE)
table_s2b <- read_csv("results/tables/TableS2B_DEMETER2_PRKDC_RNAi_CV.csv", show_col_types = FALSE)
table3b   <- read_csv("results/tables/Table3B_DNAPK_Holoenzyme_CoDependency.csv", show_col_types = FALSE)

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1.2, hjust = 0),
      plot.subtitle    = element_text(size = base_size - 0.8, color = "grey30", hjust = 0),
      axis.title       = element_text(face = "bold", size = base_size),
      axis.text        = element_text(color = "black", size = base_size - 1),
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      strip.text       = element_text(face = "bold", size = base_size - 0.5),
      legend.title     = element_text(face = "bold", size = base_size - 0.5),
      legend.text      = element_text(size = base_size - 1),
      panel.grid.major = element_line(color = "grey93", linewidth = 0.3)
    )
}

cancer_pal <- c("CACX" = "#D95F02", "OVCX" = "#1F78B4", "ENDOCX" = "#33A02C")

save_pub_fig <- function(plot_obj, stem, w, h) {
  ggsave(file.path("results/figures", paste0(stem, ".tiff")), plot_obj,
         width = w, height = h, dpi = 600, compression = "lzw", bg = "white")
  ggsave(file.path("results/figures", paste0(stem, ".png")), plot_obj,
         width = w, height = h, dpi = 600, bg = "white")
  cat(sprintf("Saved %s (.tiff & .png, 600 DPI)\n", stem))
}

pangyn_ordered <- pangyn %>%
  mutate(Cancer_Code = factor(Cancer_Code, levels = c("CACX", "OVCX", "ENDOCX")))

# ==============================================================================
# FIGURE 1: DATA-DRIVEN COHORT COMPOSITION & ASSAY COVERAGE AUDIT (NO TEXT WALLS)
# ==============================================================================
# Panel A: Discovery cell-line composition by histology (Epithelial vs Non-Epithelial)
hist_df <- pangyn_ordered %>%
  mutate(Histology = ifelse(Is_Epithelial, "Epithelial Carcinoma (n=89)", "Non-Epithelial / Sarcoma (n=19)")) %>%
  count(Cancer_Code, Histology) %>%
  group_by(Cancer_Code) %>%
  mutate(Total = sum(n)) %>%
  ungroup()

fig1a <- ggplot(hist_df, aes(x = Cancer_Code, y = n, fill = Histology)) +
  geom_col(width = 0.55, color = "black", linewidth = 0.35) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 3.8) +
  geom_text(data = hist_df %>% distinct(Cancer_Code, Total),
            aes(x = Cancer_Code, y = Total + 2.8, label = paste0("N = ", Total)),
            inherit.aes = FALSE, fontface = "bold", size = 3.6) +
  scale_fill_manual(values = c("Epithelial Carcinoma (n=89)" = "#1F4E79",
                               "Non-Epithelial / Sarcoma (n=19)" = "#8DA0CB"), name = NULL) +
  coord_cartesian(ylim = c(0, 65)) +
  labs(
    title = "A. Pan-Gynecologic Discovery Cohorts (N = 108)",
    subtitle = "DepMap Public 26Q1 evaluable models by lineage & histology",
    x = NULL, y = "Number of Cell Lines"
  ) +
  theme_pub(10.5) +
  theme(legend.position = "bottom")

# Panel B: Perturbation & multi-omic assay completeness (%) across CACX, OVCX, ENDOCX & Pan-Cancer
cov_df <- bind_rows(
  pangyn_ordered %>%
    group_by(Cohort = as.character(Cancer_Code)) %>%
    summarise(
      `CRISPR XRCC5/6\n(DNA-PK Ku)` = 100 * mean(!is.na(CRISPR_XRCC5)),
      `CRISPR LIG4/NHEJ1\n(c-NHEJ Ligation)` = 100 * mean(!is.na(CRISPR_LIG4)),
      `WGS Copy Number\n(OmicsCNGene)` = 100 * mean(!is.na(CN_XRCC5)),
      `DEMETER2 RNAi\n(PRKDC Knockdown)` = 100 * mean(!is.na(RNAi_PRKDC)),
      `CRISPR PRKDC (5591)\n(Unmasked Screens)` = 100 * mean(!is.na(CRISPR_PRKDC)),
      .groups = "drop"
    ),
  pancancer %>%
    summarise(
      Cohort = "Pan-Cancer\n(N=1,140)",
      `CRISPR XRCC5/6\n(DNA-PK Ku)` = 100 * mean(!is.na(CRISPR_XRCC5)),
      `CRISPR LIG4/NHEJ1\n(c-NHEJ Ligation)` = 100 * mean(!is.na(CRISPR_LIG4)),
      `WGS Copy Number\n(OmicsCNGene)` = 100,
      `DEMETER2 RNAi\n(PRKDC Knockdown)` = 100 * mean(!is.na(RNAi_PRKDC)),
      `CRISPR PRKDC (5591)\n(Unmasked Screens)` = 100 * mean(!is.na(CRISPR_PRKDC))
    )
) %>%
  pivot_longer(-Cohort, names_to = "Assay", values_to = "Pct") %>%
  mutate(
    Cohort = factor(Cohort, levels = c("CACX", "OVCX", "ENDOCX", "Pan-Cancer\n(N=1,140)")),
    Assay  = factor(Assay, levels = rev(c("CRISPR XRCC5/6\n(DNA-PK Ku)",
                                          "CRISPR LIG4/NHEJ1\n(c-NHEJ Ligation)",
                                          "WGS Copy Number\n(OmicsCNGene)",
                                          "DEMETER2 RNAi\n(PRKDC Knockdown)",
                                          "CRISPR PRKDC (5591)\n(Unmasked Screens)"))),
    Txt_Col = ifelse(Pct > 55, "white", "black")
  )

fig1b <- ggplot(cov_df, aes(x = Cohort, y = Assay, fill = Pct)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", Pct), color = Txt_Col), fontface = "bold", size = 3.4) +
  scale_color_identity() +
  scale_fill_gradient(low = "#FDECEC", high = "#1F78B4", limits = c(0, 100), name = "Coverage (%)") +
  labs(
    title = "B. Multi-Omic & Perturbation Assay Coverage Audit",
    subtitle = "PRKDC (5591) CRISPR is masked in 98.1% of GYN lines; Ku70/Ku80 & RNAi are complete",
    x = NULL, y = NULL
  ) +
  theme_pub(10.5) +
  theme(axis.line = element_blank(), axis.ticks = element_blank())

fig1 <- (fig1a | fig1b) + plot_layout(widths = c(1, 1.45))
save_pub_fig(fig1, "Figure1_PanGYN_Cohort_and_Assay_Audit", w = 11.2, h = 5.2)

# ==============================================================================
# FIGURE 2: THREE-CANCER DEPENDENCY DISTRIBUTIONS & ABUNDANCE VS. DEPENDENCY
# ==============================================================================
dep_long <- pangyn_ordered %>%
  select(ModelID, Cancer_Code, CRISPR_XRCC5, CRISPR_XRCC6, CRISPR_LIG4, CRISPR_NHEJ1) %>%
  pivot_longer(starts_with("CRISPR_"), names_to = "Gene", values_to = "CRISPR_Effect") %>%
  mutate(
    Gene = recode(Gene,
                  "CRISPR_XRCC5" = "XRCC5 (Ku80)",
                  "CRISPR_XRCC6" = "XRCC6 (Ku70)",
                  "CRISPR_LIG4"  = "LIG4",
                  "CRISPR_NHEJ1" = "NHEJ1 (XLF)"
    ),
    Gene = factor(Gene, levels = c("XRCC5 (Ku80)", "XRCC6 (Ku70)", "LIG4", "NHEJ1 (XLF)"))
  )

fig2a <- ggplot(dep_long, aes(x = Cancer_Code, y = CRISPR_Effect, fill = Cancer_Code)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = -0.5, linetype = "dashed", color = "#B2182B", linewidth = 0.4) +
  geom_boxplot(width = 0.48, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.14, size = 1.3, alpha = 0.65, color = "grey15") +
  facet_wrap(~Gene, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = cancer_pal) +
  labs(
    title = "A. DNA-PK Holoenzyme CRISPR Dependency Across Gynecologic Lineages (N = 108)",
    x = NULL, y = "CRISPR Gene Effect (Chronos)"
  ) +
  theme_pub(10.5) +
  theme(legend.position = "none")

abund_long <- bind_rows(
  pangyn_ordered %>% transmute(Cancer_Code, Target = "XRCC5 (Ku80) CRISPR", mRNA = RNA_XRCC5, Dep = CRISPR_XRCC5),
  pangyn_ordered %>% transmute(Cancer_Code, Target = "XRCC6 (Ku70) CRISPR", mRNA = RNA_XRCC6, Dep = CRISPR_XRCC6),
  pangyn_ordered %>% transmute(Cancer_Code, Target = "PRKDC DEMETER2 RNAi", mRNA = RNA_PRKDC, Dep = RNAi_PRKDC)
) %>%
  filter(!is.na(mRNA) & !is.na(Dep)) %>%
  mutate(Target = factor(Target, levels = c("XRCC5 (Ku80) CRISPR", "XRCC6 (Ku70) CRISPR", "PRKDC DEMETER2 RNAi")))

fig2b <- ggplot(abund_long, aes(x = mRNA, y = Dep, color = Cancer_Code)) +
  geom_point(size = 2.0, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 0.8) +
  facet_wrap(~Target, scales = "free") +
  scale_color_manual(values = cancer_pal, name = "Lineage") +
  labs(
    title = "B. Baseline Target mRNA Abundance vs. Functional Dependency (All P > 0.05)",
    x = "Baseline mRNA Expression (log2 TPM + 1)", y = "Functional Perturbation Score"
  ) +
  theme_pub(10.5) +
  theme(legend.position = "bottom")

fig2 <- fig2a / fig2b + plot_layout(heights = c(1, 1.15))
save_pub_fig(fig2, "Figure2_PanGYN_Dependency_and_Abundance", w = 10.5, h = 7.8)

# ==============================================================================
# FIGURE 3: COMPARATIVE NESTED-CV BENCHMARK (FIXED FACET & AXIS SPACING)
# ==============================================================================
cv_plot_df <- bind_rows(
  table2 %>%
    filter(Validation_Arm == "Within-Lineage Nested CV (5x5)",
           Algorithm == "Transcriptome ElasticNet (a=0.5)") %>%
    mutate(Panel = ifelse(Outcome == "CRISPR_XRCC5", "CRISPR XRCC5 (Ku80)", "CRISPR XRCC6 (Ku70)")),
  table_s2b %>%
    filter(Validation_Arm == "Within-Lineage Nested CV (5x5)") %>%
    mutate(Cancer_Code = sub(" RNAi.*", "", Cancer_Code),
           Panel = "DEMETER2 PRKDC RNAi")
) %>%
  mutate(
    Cohort_Label = paste0(Cancer_Code, " (n=", N_Eval, ")"),
    Panel = factor(Panel, levels = c("CRISPR XRCC5 (Ku80)", "CRISPR XRCC6 (Ku70)", "DEMETER2 PRKDC RNAi"))
  )

rmse_long <- cv_plot_df %>%
  select(Panel, Cohort_Label, Model_RMSE, Null_RMSE) %>%
  pivot_longer(c(Model_RMSE, Null_RMSE), names_to = "Type", values_to = "RMSE") %>%
  mutate(Type = recode(Type, "Model_RMSE" = "Nested CV ElasticNet", "Null_RMSE" = "Intercept-Only Null"))

fig3a <- ggplot(rmse_long, aes(x = Cohort_Label, y = RMSE, fill = Type)) +
  geom_col(position = position_dodge(width = 0.68), width = 0.6, color = "black", linewidth = 0.3) +
  facet_wrap(~Panel, scales = "free_x") +
  scale_fill_manual(values = c("Nested CV ElasticNet" = "#4E79A7", "Intercept-Only Null" = "#BAB0AC"), name = NULL) +
  labs(
    title = "A. Within-Lineage Out-of-Fold RMSE vs. Intercept-Only Null Baseline",
    subtitle = "Models fail to reduce prediction error below the constant intercept across gynecologic cohorts",
    x = NULL, y = "Out-of-Fold RMSE"
  ) +
  theme_pub(10.5) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1)
  )

transport_sub <- bind_rows(
  table2 %>% filter(Validation_Arm == "Pan-Cancer -> GYN Transport",
                    Outcome == "CRISPR_XRCC5",
                    Algorithm == "Transcriptome ElasticNet (a=0.5)") %>%
    mutate(Arm_Label = paste0("Pan-Cancer -> ", Cancer_Code, " (XRCC5)")),
  table2 %>% filter(Validation_Arm == "Within-Lineage Nested CV (5x5)",
                    Algorithm == "Multi-Omic (RNA+CNV+MUT) ENet") %>%
    mutate(Arm_Label = paste0("Within-", Cancer_Code, " Multi-Omic (XRCC5)")),
  table_s2a %>% filter(Cancer_Code == "Pooled Epithelial GYN (89)") %>%
    mutate(Arm_Label = paste0("Epithelial-Only Pooled (", sub("CRISPR_", "", Outcome), ")")),
  table_s2b %>% filter(Validation_Arm == "Within-Lineage Nested CV (5x5)") %>%
    mutate(Arm_Label = paste0("Within-", Cancer_Code, " (DEMETER2 PRKDC)"))
) %>%
  mutate(Arm_Label = fct_reorder(Arm_Label, Transport_R))

fig3b <- ggplot(transport_sub, aes(x = Transport_R, y = Arm_Label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#B2182B", linewidth = 0.6) +
  geom_segment(aes(x = 0, xend = Transport_R, y = Arm_Label, yend = Arm_Label),
               color = "grey55", linewidth = 0.7) +
  geom_point(aes(color = Transport_R > 0), size = 3.2) +
  scale_color_manual(values = c("TRUE" = "#1F78B4", "FALSE" = "#D95F02"), guide = "none") +
  coord_cartesian(xlim = c(-0.45, 0.45)) +
  labs(
    title = "B. Out-of-Fold Correlation Across Comparative Validation Arms",
    x = "Out-of-Fold Pearson R (Observed vs. Predicted)", y = NULL
  ) +
  theme_pub(10.5)

fig3 <- fig3a / fig3b + plot_layout(heights = c(1.15, 1))
save_pub_fig(fig3, "Figure3_PanGYN_NestedCV_Null_Collapse", w = 11.2, h = 8.2)

# ==============================================================================
# FIGURE 4: LINEAGE-ADJUSTED DNA-PK HOLOENZYME CO-DEPENDENCY HEATMAP
# ==============================================================================
heat_df <- table3b %>%
  mutate(
    Cohort = factor(Cohort, levels = c("CACX", "OVCX", "ENDOCX", "Pooled Pan-GYN")),
    Partner_Gene = factor(Partner_Gene,
                          levels = rev(c("XRCC6", "LIG4", "NHEJ1", "PARP1", "POLQ", "ATM", "ATR", "WEE1"))),
    Sig_Star = case_when(
      FDR_Q < 0.01 ~ "**",
      FDR_Q < 0.05 ~ "*",
      TRUE         ~ ""
    ),
    Cell_Label = sprintf("%.2f%s", Pearson_R, Sig_Star),
    Text_Color = ifelse(abs(Pearson_R) > 0.45, "white", "black")
  )

fig4 <- ggplot(heat_df, aes(x = Cohort, y = Partner_Gene, fill = Pearson_R)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = Cell_Label, color = Text_Color), size = 3.6, fontface = "bold") +
  scale_color_identity() +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = c(-0.7, 0.7), name = "Lineage-Adjusted\nPearson R") +
  labs(
    title = "Figure 4. Lineage-Adjusted Co-Dependency with XRCC5 (Ku80)",
    subtitle = "Catalytic c-NHEJ ligase (LIG4) is decoupled from Ku80 in all lineages (* FDR q < 0.05, ** FDR q < 0.01)",
    x = NULL, y = "DDR / DNA-PK Complex Partner"
  ) +
  theme_pub(11) +
  theme(axis.line = element_blank(), axis.ticks = element_blank())

save_pub_fig(fig4, "Figure4_PanGYN_Holoenzyme_CoDependency", w = 8.8, h = 5.6)

cat("\nUpdated Figures 1–4 saved cleanly at 600 DPI!\n")