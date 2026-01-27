library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(lme4)
#Need to run script "stats (ICC)" .... first
# 1. Read & pivot
df <- read_csv("raw_no_ctrls_with_modz.csv", show_col_types = FALSE)

df_long <- df %>%
  pivot_longer(
    cols      = starts_with("donor_"),
    names_to  = "donor",
    values_to = "value"
  ) %>%
  mutate(value = value / 100)

# 2. Compute modified‐z
df_modz <- df_long %>%
  # drop stray “...6” if present
  select(-matches("^\\.\\.\\.[0-9]+$")) %>%
  group_by(cytokine) %>%
  mutate(
    med   = median(value,            na.rm = TRUE),
    mad   = median(abs(value - med), na.rm = TRUE),
    z_mod = 0.6745 * (value - med) / mad
  ) %>%
  ungroup()

# 3. LRT for donor effect (full model first!)
donor_lrt <- df_modz %>%
  group_by(cytokine) %>%
  nest() %>%
  mutate(
    null    = map(data, ~ lm(z_mod ~ treatment, data = .x)),
    full    = map(data, ~ lmer(z_mod ~ treatment + (1 | donor),
                               data = .x, REML = FALSE)),
    p_donor = map2_dbl(full, null,
                       ~ {
                         an <- anova(.x, .y)
                         an[2, "Pr(>Chisq)"]
                       }
    )
  ) %>%
  select(cytokine, p_donor) %>%
  ungroup()

print(donor_lrt)

library(dplyr)
library(knitr)
# (optionally) library(kableExtra)  # for more styling

# 1. Assuming you already have:
#    icc_table  = tibble(cytokine, donor_var, resid_var, ICC)
#    donor_lrt  = tibble(cytokine, p_donor)

# 2. Join them
summary_df <- icc_table %>%
  left_join(donor_lrt, by = "cytokine") %>%
  select(cytokine, ICC, p_donor)

# (1) Starting from your joined summary_df:
summary_table <- summary_df %>%
  # format p‐values and recode cytokine names
  mutate(
    `p-value` = case_when(
      p_donor < 0.001 ~ "<0.001",
      TRUE            ~ sprintf("%.3f", p_donor)
    ),
    Cytokine = dplyr::recode(cytokine,
                      vegf    = "VEGF",
                      il.8    = "IL-8",
                      kc.like = "KC-like",
                      il.10   = "IL-10",
                      ccl2    = "CCL2",
                      tnf.a   = "TNF-α",
                      tgf.b   = "TGF-β"
    )
  ) %>%
  select(Cytokine, ICC, `p-value`)

# (2) Print the nicely‐formatted kable
knitr::kable(
  summary_table,
  col.names = c("Cytokine", "ICC", "*p-value*"),
  caption   = "Donor ICC and LRT p-values by Cytokine",
  digits    = c(NA, 3, NA)
)

# (3) Save to CSV
write_csv(summary_table, "donor_icc_pvalues.csv")
