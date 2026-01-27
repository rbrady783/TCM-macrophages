library(tidyverse)
library(lme4)
library(lmerTest)

# Read data
data <- read_csv("exo_count_vegf.csv")

# Convert to long format
data_long <- data %>%
  pivot_longer(
    cols = starts_with("dog_"),
    names_to = "donor",
    values_to = "VEGF"
  ) %>%
  mutate(
    donor = factor(donor),
    cell_line = factor(cell_line)
  )

# Look at the data structure
head(data_long)
summary(data_long)

# OPTION 1: Mixed-effects model (BEST - accounts for donor variation)
model1 <- lmer(VEGF ~ exo_count + (1|donor), data = data_long)
summary(model1)
anova(model1)

# OPTION 2: Z-score within donor, then analyze
data_long_z <- data_long %>%
  group_by(donor) %>%
  mutate(VEGF_z = scale(VEGF)[,1]) %>%
  ungroup()

model2 <- lm(VEGF_z ~ exo_count, data = data_long_z)
summary(model2)

# OPTION 3: Average across donors per cell line (simplest)
data_avg <- data %>%
  mutate(mean_VEGF = rowMeans(select(., starts_with("dog_")), na.rm = TRUE))

model3 <- lm(mean_VEGF ~ exo_count, data = data_avg)
summary(model3)

# Compare all three
cat("\n=== COMPARISON ===\n")
cat("Option 1 (Mixed model):\n")
print(anova(model1))
cat("\nOption 2 (Z-scored within donor):\n")
print(summary(model2)$coefficients)
cat("\nOption 3 (Averaged):\n")
print(summary(model3)$coefficients)

# For whole TCM (6 donors, z-scored)
confint(model2, "exo_count")

# For exosome-only (3 donors, z-scored) 
confint(model2, "exo_count")  # Run this on the exosome-only analysis
