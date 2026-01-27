# Load the libraries you need  
library(readr)        # for read_csv()  
library(dplyr)        # for data‐manipulation  
library(tidyr)        # for pivoting  
library(purrr)        # for map()  
library(lmerTest)     # for lmer() with p.values  
library(broom.mixed)  # for tidy()

# 1. Read your data (percent-of-median values in “3_donors_all.csv”)
df <- read_csv("raw_no_ctrls_with_modz.csv", show_col_types = FALSE)

# 2. Pivot to long form and convert to numeric fraction
df_long <- df %>%
  pivot_longer(
    cols      = starts_with("donor_"),
    names_to  = "donor",
    values_to = "value"
  ) %>%
  mutate(value = as.numeric(value) / 100)

# 3. Compute modified-z (median-based) within each cytokine
df_modz <- df_long %>%
  group_by(cytokine) %>%
  mutate(
    med   = median(value,                na.rm = TRUE),
    mad   = median(abs(value - med),     na.rm = TRUE),
    z_mod = 0.6745 * (value - med) / mad
  ) %>%
  ungroup()

# 4. Fit one mixed model per cytokine, extract donor/residual variance & ICC, plus treatment effects
results <- df_modz %>%
  group_by(cytokine) %>%
  nest() %>%
  mutate(
    model     = map(data, ~ lmer(z_mod ~ treatment + (1|donor), data = .x)),
    vc_df     = map(model, ~ as.data.frame(VarCorr(.x))),
    donor_var = map_dbl(vc_df, ~ filter(.x, grp=="donor")    %>% pull(vcov)),
    resid_var = map_dbl(vc_df, ~ filter(.x, grp=="Residual") %>% pull(vcov)),
    ICC       = donor_var / (donor_var + resid_var),
    fixed_eff = map(model, ~ tidy(.x) %>% 
                      filter(effect=="fixed", term != "(Intercept)"))
  ) %>%
  ungroup()

# 5. Extract a clean ICC table
icc_table <- results %>%
  select(cytokine, donor_var, resid_var, ICC)

print(icc_table)

# 6. Extract a tidy table of treatment fixed-effects (with estimates, SE, t-stat & p-value)
treatment_effects <- results %>%
  select(cytokine, fixed_eff) %>%
  unnest(fixed_eff) %>%
  select(cytokine, term, estimate, std.error, statistic, p.value)

print(treatment_effects, n=175)
