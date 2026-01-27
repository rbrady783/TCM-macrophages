# Load libraries
library(tidyverse)
library(lme4)
library(lmerTest)

# 1. READ AND PREPARE DATA
df <- read_csv("raw_no_ctrls.csv", show_col_types = FALSE)

# Pivot to long format
df_long <- df %>%
  pivot_longer(
    cols = starts_with("donor_"),
    names_to = "donor",
    values_to = "value"
  ) %>%
  mutate(
    value = as.numeric(value),
    donor = factor(donor)
  )

# 2. COMPUTE MODIFIED Z-SCORES (using your exact formula)
df_modz <- df_long %>%
  group_by(cytokine) %>%
  mutate(
    med_val = median(value, na.rm = TRUE),
    mad_val = median(abs(value - med_val), na.rm = TRUE),
    z_mod = case_when(
      mad_val == 0 ~ 0,
      TRUE ~ 0.6745 * (value - med_val) / mad_val
    )
  ) %>%
  ungroup()

# 3. CALCULATE MEAN Z-SCORE PER TREATMENT, THEN MEDIAN & RANGE
# This matches what you show in your Word table
descriptive_stats <- df_modz %>%
  group_by(cytokine, treatment) %>%
  summarise(mean_z = mean(z_mod, na.rm = TRUE), .groups = "drop") %>%
  group_by(cytokine) %>%
  summarise(
    median_z = median(mean_z, na.rm = TRUE),
    min_z = min(mean_z, na.rm = TRUE),
    max_z = max(mean_z, na.rm = TRUE),
    range_z = max_z - min_z,
    .groups = "drop"
  )

print("Descriptive statistics (median of mean z-scores per treatment):")
print(descriptive_stats)

# 4. CALCULATE ICCs AND TEST DONOR EFFECTS
print("\nCalculating ICCs and chi-square tests...")

icc_results <- df_modz %>%
  group_by(cytokine) %>%
  nest() %>%
  mutate(
    # Fit models
    model_full = map(data, ~ lmer(z_mod ~ treatment + (1|donor), 
                                  data = .x, REML = FALSE)),
    model_null = map(data, ~ lm(z_mod ~ treatment, data = .x)),
    
    # Extract variance components for ICC
    vc_df = map(model_full, ~ as.data.frame(VarCorr(.x))),
    donor_var = map_dbl(vc_df, ~ filter(.x, grp == "donor") %>% pull(vcov)),
    resid_var = map_dbl(vc_df, ~ filter(.x, grp == "Residual") %>% pull(vcov)),
    ICC = donor_var / (donor_var + resid_var),
    
    # Likelihood ratio test for donor effect
    lrt = map2(model_full, model_null, ~ anova(.x, .y)),
    p_donor = map_dbl(lrt, ~ .x[2, "Pr(>Chisq)"]),
    chi_sq = map_dbl(lrt, ~ .x[2, "Chisq"])
  ) %>%
  ungroup()

# 5. ADD BH ADJUSTMENT AND JOIN DESCRIPTIVE STATS
icc_results <- icc_results %>%
  mutate(p_donor_adj = p.adjust(p_donor, method = "BH")) %>%
  left_join(descriptive_stats, by = "cytokine")

# 6. CREATE SUMMARY TABLE (sorted by median, no redundant columns)
summary_table <- icc_results %>%
  select(cytokine, median_z, range_z, ICC, chi_sq, p_donor_adj) %>%
  arrange(desc(median_z)) %>%
  mutate(
    Cytokine = dplyr::recode(cytokine,
                             "vegf" = "VEGF",
                             "il.8" = "IL-8",
                             "kc.like" = "KC-like",
                             "il.10" = "IL-10",
                             "ccl2" = "CCL2",
                             "tnf.a" = "TNF-α",
                             "tgf.b" = "TGF-β"),
    Median = sprintf("%.3f", median_z),
    Range = sprintf("%.3f", range_z),
    ICC = sprintf("%.3f", ICC),
    `χ²` = sprintf("%.2f", chi_sq),
    `Adj. p-value` = case_when(
      p_donor_adj < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_donor_adj)
    )
  ) %>%
  select(Cytokine, Median, Range, ICC, `χ²`, `Adj. p-value`)

print("\nFinal Summary Table:")
print(summary_table)

# Save results
write_csv(summary_table, "donor_icc_results_complete.csv")
