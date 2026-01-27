library(tidyverse)

# Read data
data <- read_csv("rna_protein.csv")

# Clean column names
colnames(data) <- c("cell_line", "CCL3_protein", "CCL3_RNA")

# Summary statistics
cat("=== CCL3 PROTEIN SUMMARY ===\n")
summary(data$CCL3_protein)

# Identify groups
data <- data %>%
  mutate(
    group = case_when(
      cell_line %in% c("DH82", "Nike") ~ "High (HS)",
      CCL3_protein > 0 ~ "Low detectable",
      TRUE ~ "Undetectable"
    )
  )

# Group summaries
group_summary <- data %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean = mean(CCL3_protein),
    sd = sd(CCL3_protein),
    min = min(CCL3_protein),
    max = max(CCL3_protein)
  )

print(group_summary)

# DH82 and Nike specific values
dh82_nike <- data %>% 
  filter(cell_line %in% c("DH82", "Nike"))
print(dh82_nike)

# Detection limit
cat("\nDetection limit (minimum non-zero):", min(data$CCL3_protein[data$CCL3_protein > 0]), "pg/mL\n")
