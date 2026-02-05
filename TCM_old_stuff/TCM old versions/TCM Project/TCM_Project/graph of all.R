library(dplyr)
library(tidyr)

# Read the CSV file (using check.names = FALSE preserves original column names)
df <- read.csv("tumortype.csv", check.names = FALSE)

# Rename columns for easier reference.
# Here, we rename "Cell line" to "Cell_line" and "Tumor.type" to "TumorType".
df <- df %>%
  rename(Cell_line = `Cell line`,
         TumorType = `Tumor.type`)

# Pivot the Dog columns from wide to long format.
# This creates two new columns: "Dog" (the replicate identifier) and "Value" (the measurement)
df_tidy <- df %>%
  pivot_longer(
    cols = starts_with("Dog"),   # this grabs "Dog 1", "Dog 2", "Dog 3"
    names_to = "Dog",
    values_to = "Value"
  )

# Inspect the tidy data:
head(df_tidy)

summary_df <- df_tidy %>%
  group_by(Cytokine, TumorType) %>%
  summarise(
    mean_value = mean(Value, na.rm = TRUE),
    se = sd(Value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

cytokine_order <- df_tidy %>%
  group_by(Cytokine) %>%
  summarise(overall_mean = mean(Value, na.rm = TRUE)) %>%
  arrange(desc(overall_mean)) %>%
  pull(Cytokine)

# Reorder the Cytokine factor in the summary data frame:
summary_df <- summary_df %>%
  mutate(Cytokine = factor(Cytokine, levels = cytokine_order))

library(ggplot2)

ggplot(summary_df, aes(x = Cytokine, y = mean_value, fill = TumorType)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_value - se, ymax = mean_value + se), width = 0.2) +
  labs(
    x = "Cytokine",
    y = "Mean Value ± SE",
    title = "Cytokine Means by Tumor Type"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
