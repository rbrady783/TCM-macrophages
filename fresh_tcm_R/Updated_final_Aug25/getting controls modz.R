library(readr)
library(dplyr)
library(tidyr)
library(tibble)

# Read original data (not z-scores)
df <- read_csv("raw_no_ctrls.csv", show_col_types = FALSE)

# Pivot to long format
df_long <- df %>%
  pivot_longer(
    cols = starts_with("donor_"),
    names_to = "donor",
    values_to = "value"
  ) %>%
  mutate(value = as.numeric(value))

# Calculate median & MAD from the 25 treatments (excluding controls)
cyto_stats <- df_long %>%
  group_by(cytokine) %>%
  summarise(
    med_val = median(value, na.rm = TRUE),
    mad_val = median(abs(value - med_val), na.rm = TRUE),
    .groups = "drop"
  )

# Your control data
df_control <- tibble(
  cytokine  = c("vegf", "il.8", "kc.like", "il.10", "ccl2", "tnf.a", "tgf.b"),
  treatment = "Ctrl",
  donor_1   = c(4254.387207, 8502, 607.33, 56.87, 18073.48688, 6.9, 18.8012),
  donor_2   = c(1547.878878, 22561, 774.74, 105.96, 222748.8108, 6.64, 37.2757),
  donor_3   = c(1749.378786, 29578, 939.54, 82.65, 116809.5937, 6.87, 137.2803)
)

# Pivot controls to long form
df_control_long <- df_control %>%
  pivot_longer(
    cols = starts_with("donor_"),
    names_to = "donor",
    values_to = "value"
  )

# Compute modified z-scores for controls using treatment-derived stats
df_control_modz <- df_control_long %>%
  left_join(cyto_stats, by = "cytokine") %>%
  mutate(
    z_modified_control = case_when(
      mad_val == 0 ~ 0,
      TRUE ~ 0.6745 * (value - med_val) / mad_val
    )
  )

# Check the results
print(df_control_modz %>% 
        select(cytokine, donor, value, med_val, mad_val, z_modified_control))

# Pivot back to wide format
df_control_modz_wide <- df_control_modz %>%
  select(cytokine, treatment, donor, z_modified_control) %>%
  pivot_wider(
    names_from = donor,
    values_from = z_modified_control,
    names_glue = "{donor}_modz"
  )

print(df_control_modz_wide)
write_csv(df_control_modz_wide, "ctrls_modz.csv")


      