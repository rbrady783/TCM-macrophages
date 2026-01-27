# Statistical Analysis: MVB12A High vs Low Expression
# Load required libraries
library(tidyverse)
library(broom)

# Read the data
data <- read.csv("mvb12A_with_CIN.csv")

# Clean up column names and inspect data
colnames(data) <- c("Cell_Line", "Dog1", "Dog2", "Dog3", "MVB12A")
head(data)
str(data)

# Convert MVB12A to factor
data$MVB12A <- as.factor(data$MVB12A)
summary(data)


# ============================================================================
# APPROACH 2: Average across dogs (single t-test)
# ============================================================================

print("\n=== APPROACH 2: Average response across dogs ===")

# Calculate mean response for each cell line
data$Mean_Response <- rowMeans(data[,c("Dog1", "Dog2", "Dog3")])

# Descriptive statistics
desc_stats <- data %>%
  group_by(MVB12A) %>%
  summarise(
    n = n(),
    mean = round(mean(Mean_Response), 1),
    sd = round(sd(Mean_Response), 1),
    se = round(sd(Mean_Response)/sqrt(n()), 1),
    .groups = 'drop'
  )
print(desc_stats)

# Welch's t-test (unequal variances)
avg_test <- t.test(Mean_Response ~ MVB12A, data = data)
print(avg_test)

# Effect size (Cohen's d)
high_vals <- data$Mean_Response[data$MVB12A == "high"]
low_vals <- data$Mean_Response[data$MVB12A == "low"]
pooled_sd <- sqrt(((length(high_vals)-1)*sd(high_vals)^2 + (length(low_vals)-1)*sd(low_vals)^2) / 
                    (length(high_vals) + length(low_vals) - 2))
cohens_d <- (mean(high_vals) - mean(low_vals)) / pooled_sd
cat("Cohen's d (effect size):", round(cohens_d, 2), "\n")

##normal data
# Check residuals from your t-test
residuals <- c(high_vals - mean(high_vals), low_vals - mean(low_vals))
shapiro.test(residuals)  # Test normality
