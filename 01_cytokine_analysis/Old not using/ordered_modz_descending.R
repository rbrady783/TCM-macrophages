#Getting mod z scores in descending order

library(dplyr)
library(readr)

# 1. Read your already‐wide CSV
df <- read_csv("modz_with_mean.csv", show_col_types = FALSE)

# 2. Reorder each cytokine’s treatments by descending mean
df_ordered <- df %>%
  group_by(cytokine) %>%
  arrange(desc(mean), .by_group = TRUE) %>%
  ungroup()

# 3. (Optional) drop the mean column if you just want the replicates
df_ordered %>% 
  select(cytokine, treatment, donor_1_modz, donor_2_modz, donor_3_modz) %>%
  print(n = Inf)

# 4. (Optional) write it back out
write_csv(df_ordered, "modz_replicates_ordered_by_cytokine.csv")
