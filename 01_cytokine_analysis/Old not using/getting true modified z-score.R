library(readr)
library(dplyr)
library(tidyr)

# 1. Read in your raw data
df <- read_csv("C:/Users/brady/OneDrive/Desktop/fresh_tcm_R/Updated_final_Aug25/raw_no_ctrls.csv")

# 2. Pivot to long form so you have one row per cytokine × treatment × donor
df_long <- df %>%
  pivot_longer(
    cols      = starts_with("donor_"),
    names_to  = "donor",
    values_to = "value"
  ) %>%
  mutate(value = as.numeric(value))

# 3. Compute z-scores _within each cytokine_
df_z <- df_long %>%
  group_by(cytokine) %>%
  mutate(
    # standard z‐score: (x - mean) / sd
    mean_val    = mean(value, na.rm = TRUE),
    sd_val      = sd(value, na.rm = TRUE),
    z_standard  = (value - mean_val) / sd_val,
    
    # modified z‐score: 0.6745 * (x - median) / MAD
    med_val     = median(value, na.rm = TRUE),
    mad_val     = median(abs(value - med_val), na.rm = TRUE),
    z_modified  = 0.6745 * (value - med_val) / mad_val
  ) %>%
  ungroup() %>%
  select(cytokine, treatment, donor, value, z_standard, z_modified)

# 4. Take a peek
print(df_z %>% slice_head(n = 10))


# 4. Pivot back to wide form so each donor has its own modified-z column
df_modz_wide <- df_z %>%
  select(cytokine, treatment, donor, z_modified) %>%
  pivot_wider(
    names_from  = donor,
    values_from = z_modified,
    names_glue  = "{donor}_modz"
  )

# 5. Inspect the first few rows
print(df_modz_wide %>% slice_head(n = 10))

# (Optional) write the result out to a new CSV
write_csv(df_modz_wide, "raw_no_ctrls_with_modz.csv")
