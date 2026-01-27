# Load required libraries
library(lme4)
library(lmerTest)
library(dplyr)
library(tidyr)

# Read the data from CSV file
data <- read.csv("vegf_in_tcm_vs_macs.csv")

# Reshape data from wide to long format
full_data <- data %>%
  pivot_longer(
    cols = c(`X1`, `X2`, `X3`, `X4`, `X5`, `X6`),
    names_to = "Donor",
    values_to = "Macrophage_VEGF"
  ) %>%
  rename(
    Cell_Line = cell_line,
    VEGF_TCM = vegf_cell
  ) %>%
  mutate(
    Donor = paste0("Donor_", Donor),
    VEGF_TCM_scaled = scale(VEGF_TCM)[,1]  # standardized
  )

# Fit mixed-effects model with scaled predictor
model <- lmer(Macrophage_VEGF ~ VEGF_TCM_scaled + (1|Donor), data = full_data)

# Results
summary(model)
anova(model)

# Calculate 95% CI
confint(model, method = "Wald")
