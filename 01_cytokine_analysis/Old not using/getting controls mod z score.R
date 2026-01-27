library(readr)
library(dplyr)
library(tidyr)
library(tibble)

#(ned to rn getting true modified z score first)
# 6. Recompute each cytokine’s median & MAD from your original data (excluding controls)
cyto_stats <- df_z %>%
  group_by(cytokine) %>%
  summarise(
    med_val = median(value,                na.rm = TRUE),
    mad_val = median(abs(value - med_val), na.rm = TRUE),
    .groups = "drop"
  )

# 7. Enter your control measurements in the same wide format 187.7802547	207.1426061	224.8556178

df_control <- tibble(
  cytokine  = c("vegf",    "il.8",     "kc.like", "il.10",       "ccl2",          "tnf.a", "tgf.b"),
  treatment = "Ctrl",
  donor_1   = c(4254.387207, 8502,   607.33,     56.87,        18073.48688,      6.9,    187.7803),
  donor_2   = c(1547.878878,22561,   774.74,    105.96,       222748.8108,    6.64,     207.1426),
  donor_3   = c(1749.378786,29578,   939.54,     82.65,       116809.5937,     6.87,    224.8556)
)

# 8. Pivot controls into long form
df_control_long <- df_control %>%
  pivot_longer(
    cols      = starts_with("donor_"),
    names_to  = "donor",
    values_to = "value"
  )

# 9. Join with cyto_stats and compute modified z-score for each control
df_control_modz <- df_control_long %>%
  left_join(cyto_stats, by = "cytokine") %>%
  mutate(
    z_modified_control = 0.6745 * (value - med_val) / mad_val
  )

# 10. (Optional) Pivot back to wide so each donor’s control z shows in its own column
df_control_modz_wide <- df_control_modz %>%
  select(cytokine, treatment, donor, z_modified_control) %>%
  pivot_wider(
    names_from  = donor,
    values_from = z_modified_control,
    names_glue  = "{donor}_modz"
  )

# Inspect
print(df_control_modz)
print(df_control_modz_wide)

write_csv(df_control_modz_wide, "ctrls_modz.csv")
