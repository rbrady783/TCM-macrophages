#for part II- look at impact of cell counts on cytokine secretion

library(tidyverse)    
library(lme4)         
library(lmerTest)     
library(broom.mixed)
library(performance)


# 1) Read & pivot RAW cytokine data
raw_df <- read_csv("raw_no_ctrls.csv")

# pivot

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

# 2) Read cell-counts table
cell_counts <- read_csv("cell_counts.csv") %>%
  rename(treatment = cell_line) %>%
  mutate(treatment = str_trim(treatment))

#(name mismatches)
raw_long <- raw_long %>%
  mutate(
    treatment = recode(treatment,
                       "D-17"    = "D17",
                       "CML-6M"  = "CML-6m",
                       "DEN-HSA" = "Denver"))


# 3) Join by treatment
df_cc_raw <- raw_long %>%
  left_join(cell_counts, by = "treatment")

# 4) Check to make sure no missing values since I fixed the names
df_cc_raw %>% 
  summarise(
    total_rows    = n(),
    missing_cells = sum(is.na(cell_count)),
    missing_raw   = sum(is.na(raw_val))
  ) %>% print()


# 5) Take a look
df_cc_raw %>% 
  select(cytokine, treatment, donor, raw_val, cell_count) %>% 
  slice(1:10) %>% 
  print(n = 10)

# Now can do per‐cytokine correlation

cell_count_results <- df_cc_raw %>%
  mutate(cell_millions = cell_count / 1e6) %>%
  group_by(cytokine) %>%
  nest() %>%
  mutate(
    fit = map(data, ~ lmer(raw_val ~ cell_millions + (1 | donor), data = .x)),
    tidy = map(fit, ~ tidy(.x, effects = "fixed"))
  ) %>%
  unnest(tidy) %>%
  filter(term == "cell_millions") %>%
  select(
    Cytokine  = cytokine,
    Estimate  = estimate,
    `Std. Error` = std.error,
    `p-value` = p.value)

# Print kable
kable(
  cell_count_results,
  caption = "Linear mixed‐model slope of cytokine vs. cell count (per million cells)",
  digits  = c(NA, 3, 3, 3)
)

# spot‐check assumptions(il10 and tnfa had the biggest spread)
il10_fit <- df_cc_raw %>%
  mutate(cell_millions = cell_count / 1e6) %>%
  filter(cytokine == "il.10") %>%
  lmer(raw_val ~ cell_millions + (1 | donor), data = .)

check_model(il10_fit, check = c("linearity", "homogeneity", "qq"))

