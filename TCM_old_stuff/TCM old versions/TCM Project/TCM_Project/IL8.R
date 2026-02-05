library(tidyr)
library(dplyr)
library(ggplot2)

# Step 1: Read and tidy the data.
# Read the CSV file; check.names = FALSE preserves original column names.
df <- read.csv("tumortype.csv", check.names = FALSE)

# Clean the TumorType column by trimming whitespace
df$Tumor.type <- trimws(df$Tumor.type)

# Rename columns for easier reference. (Assuming the first column is "Cell line" and tumor type is "Tumor.type".)
df <- df %>%
  rename(Cell_line = `Cell line`,
         TumorType = `Tumor.type`)

# Pivot the Dog columns into long format.
df_tidy <- df %>%
  pivot_longer(
    cols = starts_with("Dog"),  # picks up "Dog 1", "Dog 2", "Dog 3"
    names_to = "Dog",           # creates a new column for the dog identifier
    values_to = "Value"         # creates a new column for the measurement
  )

# Step 2: Filter for a single cytokine, e.g., "IL.8"
df_cytokine <- df_tidy %>%
  filter(Cytokine == "IL.8")


# Summarize the data by Cell_line and TumorType: compute the mean and standard error (SEM)
summary_df <- df_cytokine %>%
  group_by(Cell_line, TumorType) %>%
  summarise(
    mean_value = mean(Value, na.rm = TRUE),
    se = sd(Value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Order the Cell_line factor by descending mean_value
summary_df <- summary_df %>%
  arrange(desc(mean_value)) %>%
  mutate(Cell_line = factor(Cell_line, levels = unique(Cell_line)))

# Create the bar graph with error bars and overlay individual points for the three values
ggplot(summary_df, aes(x = Cell_line, y = mean_value, fill = TumorType)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_value - se, ymax = mean_value + se), width = 0.2) +
  # Overlay the individual data points using the original data for IL.8
  geom_jitter(data = df_cytokine, 
              aes(x = Cell_line, y = Value),
              width = 0.2, height = 0, shape = 21, size = 3, 
              fill = "white", color = "black") +
  scale_fill_manual(values = c(
    "Osteosarcoma"        = "blue",
    "Leukemia/Lymphoma"   = "red",
    "Melanoma"            = "gray",
    "TCC"                 = "green",
    "Mammary"             = "palegreen",
    "Thyroid carcinoma"   = "seagreen",
    "Hemangiosarcoma"     = "deepskyblue",
    "Histiocytic sarcoma" = "tomato1",
    "Soft tissue sarcoma" = "cornflowerblue"
  )) +
  labs(title = "IL-8",
    fill = "Tumor Type"
  ) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),  # Remove x-axis title
         axis.title.y = element_blank(),  # Remove y-axis title
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14))
