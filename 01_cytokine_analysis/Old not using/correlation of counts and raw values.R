library(tidyverse)

# 1) Read & pivot your RAW cytokine data
raw_df <- read_csv("raw_no_ctrls.csv")

#   cols: cytokine | treatment | donor_1 | donor_2 | donor_3

raw_long <- raw_df %>%
  pivot_longer(
    cols         = starts_with("donor_"),
    names_to     = "donor",
    names_prefix = "donor_",
    values_to    = "raw_val"
  ) %>%
  mutate(
    donor     = str_trim(donor),
    treatment = str_trim(treatment),
    cytokine  = str_trim(cytokine)
  )

# 2) Read your cell-counts table
cell_counts <- read_csv("cell_counts.csv") %>%
  rename(treatment = cell_line) %>%
  mutate(treatment = str_trim(treatment))

# 3) Join them by treatment
df_cc_raw <- raw_long %>%
  left_join(cell_counts, by = "treatment")

# 4) Sanity‐check
df_cc_raw %>% 
  summarise(
    total_rows    = n(),
    missing_cells = sum(is.na(cell_count)),
    missing_raw   = sum(is.na(raw_val))
  ) %>% print()


# 5) Peek at the result
df_cc_raw %>% 
  select(cytokine, treatment, donor, raw_val, cell_count) %>% 
  slice(1:10) %>% 
  print(n = 10)

# e.g. per‐cytokine correlation
library(tidyverse)
library(lme4)
library(broom.mixed)

# 1) Make sure you have df_cc_raw:
#    a table with columns: cytokine, treatment, donor, raw_val, cell_count
library(lmerTest)    # instead of lme4
library(broom)

cell_effects_scaled <- df_cc_raw %>%
  mutate(cell_millions = cell_count/1e6) %>%
  group_by(cytokine) %>%
  nest() %>%
  mutate(
    fit  = map(data, ~ lmer(raw_val ~ cell_millions + (1|donor), data = .x)),
    tab  = map(fit,  ~ broom::tidy(.x, effects="fixed"))
  ) %>%
  unnest(tab) %>%
  filter(term == "cell_millions") %>%
  select(cytokine, estimate, std.error, p.value)

print(cell_effects_scaled)

library(dplyr)
library(knitr)

# Prepare the table for kable
cell_table <- cell_effects_scaled %>%
  # format p‐values
  mutate(
    `p-value` = case_when(
      p.value < 0.001 ~ "<0.001",
      TRUE            ~ sprintf("%.3f", p.value)
    )
  ) %>%
  select(
    Cytokine = cytokine,
    Estimate = estimate,
    `Std. Error` = std.error,
    `p-value`
  )

# Print as a kable
knitr::kable(
  cell_table,
  caption = "Linear mixed‐model slope of cytokine vs. cell count (per million cells)",
  digits  = c(NA, 3, 3, NA)
)

write_csv(cell_table, "cell_count_correlation.csv")
