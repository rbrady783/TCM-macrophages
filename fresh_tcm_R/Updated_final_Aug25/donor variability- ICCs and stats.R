# Load libraries
library(tidyverse)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(RColorBrewer)

# 1. READ AND PREPARE DATA
# Read the ORIGINAL data (not the one with z-scores already computed)
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
    donor = factor(donor)  # Make donor a factor for mixed models
  )

# 2. COMPUTE MODIFIED Z-SCORES (once, correctly)
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

# 3. CALCULATE ICCs AND TEST DONOR EFFECTS
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
    p_donor = map_dbl(lrt, ~ .x[2, "Pr(>Chisq)"])
  ) %>%
  ungroup()

# Add this after your icc_results calculation:
icc_results <- icc_results %>%
  mutate(p_donor_adj = p.adjust(p_donor, method = "BH"))

# 4. CREATE SUMMARY TABLE
summary_table <- icc_results %>%
  select(cytokine, ICC, p_donor_adj) %>%  # Changed from p_donor to p_donor_adj
  mutate(
    Cytokine = dplyr::recode(cytokine,
                             "vegf" = "VEGF",
                             "il.8" = "IL-8",
                             "kc.like" = "KC-like",
                             "il.10" = "IL-10",
                             "ccl2" = "CCL2",
                             "tnf.a" = "TNF-α",
                             "tgf.b" = "TGF-β"),
    `p-value` = case_when(
      p_donor_adj < 0.001 ~ "<0.001",  # Changed to p_donor_adj
      p_donor_adj < 0.01 ~ sprintf("%.3f**", p_donor_adj),  # Changed
      p_donor_adj < 0.05 ~ sprintf("%.3f*", p_donor_adj),  # Changed
      TRUE ~ sprintf("%.3f", p_donor_adj)  # Changed
    )
  ) %>%
  select(Cytokine, ICC, `p-value`)

# Print results
print(summary_table)

# Save results
write_csv(summary_table, "donor_icc_results.csv")

# 5. VISUALIZATION - Inter-donor variability
# Calculate mean and SD across donors for each cytokine
donor_variability <- df_modz %>%
  group_by(cytokine, donor) %>%
  summarise(donor_mean = mean(z_mod, na.rm = TRUE), .groups = "drop") %>%
  group_by(cytokine) %>%
  summarise(
    mean_across_donors = mean(donor_mean),
    sd_across_donors = sd(donor_mean),
    .groups = "drop"
  ) %>%
  mutate(
    cytokine = case_when(
      cytokine == "vegf" ~ "VEGF",
      cytokine == "il.8" ~ "IL-8",
      cytokine == "kc.like" ~ "KC-like",
      cytokine == "il.10" ~ "IL-10",
      cytokine == "ccl2" ~ "CCL2",
      cytokine == "tnf.a" ~ "TNF-α",
      cytokine == "tgf.b" ~ "TGF-β",
      TRUE ~ cytokine
    ),
    cytokine = factor(cytokine, levels = c("VEGF", "IL-8", "IL-10", 
                                           "TNF-α", "TGF-β", "KC-like", "CCL2"))
  )

# Create plot (did not end up using plot, just reported in table)
palette_cols <- brewer.pal(n = 7, name = "Set2")

ggplot(donor_variability, aes(x = cytokine, y = mean_across_donors, fill = cytokine)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_across_donors - sd_across_donors, 
                    ymax = mean_across_donors + sd_across_donors),
                width = 0.2, linewidth = 0.8) +
  scale_fill_manual(values = palette_cols) +
  labs(
    title = "Inter-donor Variability in Cytokine Expression",
    subtitle = "Mean ± SD of donor means (n = 3 donors)",
    x = NULL,
    y = "Modified Z-score",
    caption = paste0("ICC values indicate proportion of variance due to donor.\n",
                     "Higher SD bars = more inter-donor variability")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, margin = margin(b = 10)),
    axis.text.x = element_text(angle = 25, hjust = 1),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "none"
  )

# Save plot
ggsave("donor_variability_plot.pdf", width = 8, height = 6)

