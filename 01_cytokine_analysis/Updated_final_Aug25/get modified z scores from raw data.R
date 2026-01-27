library(readr)
library(dplyr)
library(tidyr)

# 1. Read in your raw data
df <- read_csv("C:/Users/brady/OneDrive/Desktop/fresh_tcm_R/Updated_final_Aug25/raw_no_ctrls.csv")

# 2. Pivot to long form
df_long <- df %>%
  pivot_longer(
    cols = starts_with("donor_"),
    names_to = "donor",
    values_to = "value"
  ) %>%
  mutate(value = as.numeric(value))

# 3. Compute modified z-scores within each cytokine
df_z <- df_long %>%
  group_by(cytokine) %>%
  mutate(
    # Calculate median and MAD manually (as you did)
    med_val = median(value, na.rm = TRUE),
    mad_val = median(abs(value - med_val), na.rm = TRUE),
    
    # Modified z-score using the correct formula
    z_modified = case_when(
      mad_val == 0 ~ 0,  # Handle case where MAD = 0 (all values identical)
      TRUE ~ 0.6745 * (value - med_val) / mad_val
    )
  ) %>%
  ungroup() %>%
  select(cytokine, treatment, donor, value, z_modified)

# 4. Check for any extreme outliers (optional but useful)
outliers <- df_z %>%
  filter(abs(z_modified) > 3.5) %>%
  arrange(desc(abs(z_modified)))

if(nrow(outliers) > 0) {
  cat("Found", nrow(outliers), "extreme outliers (|z| > 3.5):\n")
  print(outliers)
}

# 5. Pivot back to wide form
df_modz_wide <- df_z %>%
  select(cytokine, treatment, donor, z_modified) %>%
  pivot_wider(
    names_from = donor,
    values_from = z_modified,
    names_glue = "{donor}_modz"
  )

# 6. Save the result
write_csv(df_modz_wide, "raw_no_ctrls_with_modz_v3.csv")

# Optional: Summary statistics
cat("\nSummary of modified z-scores:\n")
print(summary(df_z$z_modified))
