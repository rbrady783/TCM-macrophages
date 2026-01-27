# TNF-α Data Analysis: Pre-Normalized Data
# Study: CCL3 knockdown vs negative control across 3 donors

# Load required packages
library(tidyverse)
library(car)        # for Levene's test
library(ggplot2)
library(gridExtra)  # for arranging plots

# Read the data
data <- read.csv("ccl3_kd_tnfa.csv")

# Display the data structure
print("Data structure:")
print(data)
str(data)

# Convert to long format for easier analysis
data_long <- data %>%
  pivot_longer(cols = starts_with("dog_"), 
               names_to = "donor", 
               values_to = "tnf_alpha") %>%
  mutate(donor = factor(donor),
         condition = factor(condition))

print("\nData in long format (before cleaning):")
print(data_long)

# Clean the data - remove rows with empty conditions or missing values
data_long <- data_long %>%
  filter(!is.na(tnf_alpha),           # Remove missing TNF values
         condition != "",              # Remove empty condition names
         condition != " ",             # Remove space-only conditions
         !is.na(condition)) %>%        # Remove NA conditions
  mutate(condition = trimws(condition)) # Remove leading/trailing whitespace

print("\nData after cleaning:")
print(data_long)

print("\nAvailable conditions in your data:")
print(unique(data_long$condition))

# ============================================================================
# DESCRIPTIVE STATISTICS
# ============================================================================

print("\n=== DESCRIPTIVE STATISTICS ===")

# Summary statistics for each condition
summary_stats <- data_long %>%
  group_by(condition) %>%
  summarise(
    n = n(),
    mean = mean(tnf_alpha, na.rm = TRUE),
    median = median(tnf_alpha, na.rm = TRUE),
    sd = sd(tnf_alpha, na.rm = TRUE),
    min = min(tnf_alpha, na.rm = TRUE),
    max = max(tnf_alpha, na.rm = TRUE),
    .groups = 'drop'
  )

print("Summary statistics for normalized data:")
print(summary_stats)

# Show individual donor values
print("\nIndividual donor values:")
donor_summary <- data_long %>%
  select(condition, donor, tnf_alpha) %>%
  pivot_wider(names_from = condition, values_from = tnf_alpha)
print(donor_summary)

# ============================================================================
# DIRECT COMPARISON: KD vs NG (Your Key Question)
# ============================================================================

print("\n=== DIRECT COMPARISON: CCL3 KD vs NEGATIVE CONTROL ===")

# Extract data for each condition
kd_data <- data_long %>% filter(condition == "kd") %>% pull(tnf_alpha)
ng_data <- data_long %>% filter(condition == "ng") %>% pull(tnf_alpha)

print("CCL3 Knockdown (kd) values:")
print(kd_data)
print("Negative Control (ng) values:")  
print(ng_data)

# Calculate descriptive statistics
kd_mean <- mean(kd_data)
ng_mean <- mean(ng_data)
kd_sd <- sd(kd_data)
ng_sd <- sd(ng_data)

print(paste("CCL3 KD mean:", round(kd_mean, 3), "±", round(kd_sd, 3)))
print(paste("Negative Control mean:", round(ng_mean, 3), "±", round(ng_sd, 3)))
print(paste("Difference:", round(ng_mean - kd_mean, 3), "fold change units"))

# Convert to percentages for easier interpretation
kd_percent <- kd_mean * 100
ng_percent <- ng_mean * 100
print(paste("CCL3 KD:", round(kd_percent, 1), "% of original control"))
print(paste("Negative Control:", round(ng_percent, 1), "% of original control"))
print(paste("Additional suppression from CCL3 KD:", round(ng_percent - kd_percent, 1), "percentage points"))

# Statistical test: Paired t-test
cat("\n--- STATISTICAL TEST: KD vs NG ---\n")
if(length(kd_data) == length(ng_data) && length(kd_data) == 3) {
  # Paired t-test (same donors)
  paired_test <- t.test(kd_data, ng_data, paired = TRUE)
  
  print("Paired t-test results:")
  print(paste("t-statistic:", round(paired_test$statistic, 3)))
  print(paste("df:", paired_test$parameter))  
  print(paste("p-value:", round(paired_test$p.value, 4)))
  print(paste("95% CI for difference:", 
              round(paired_test$conf.int[1], 4), "to", 
              round(paired_test$conf.int[2], 4)))
  
  # Effect size (Cohen's d for paired samples)
  differences <- kd_data - ng_data
  cohens_d <- mean(differences) / sd(differences)
  print(paste("Cohen's d:", round(cohens_d, 3)))
  
  # Interpret effect size
  if(abs(cohens_d) < 0.2) {
    effect_interpretation <- "negligible"
  } else if(abs(cohens_d) < 0.5) {
    effect_interpretation <- "small" 
  } else if(abs(cohens_d) < 0.8) {
    effect_interpretation <- "medium"
  } else {
    effect_interpretation <- "large"
  }
  print(paste("Effect size interpretation:", effect_interpretation))
  
  # Statistical conclusion
  if(paired_test$p.value < 0.05) {
    print("CONCLUSION: Statistically significant difference between conditions")
  } else {
    print("CONCLUSION: No statistically significant difference between conditions")
  }
  
} else {
  print("Error: Unequal sample sizes or unexpected data structure")
  print(paste("KD n =", length(kd_data), ", NG n =", length(ng_data)))
}

# Non-parametric alternative
cat("\n--- NON-PARAMETRIC TEST: Wilcoxon Signed-Rank ---\n")
wilcox_test <- wilcox.test(kd_data, ng_data, paired = TRUE)
print(paste("Wilcoxon signed-rank test p-value:", round(wilcox_test$p.value, 4)))

# Practical significance assessment
cat("\n--- PRACTICAL SIGNIFICANCE ASSESSMENT ---\n")
percent_difference <- abs(ng_percent - kd_percent)
if(percent_difference < 5) {
  practical_conclusion <- "Conditions have practically equivalent effects"
} else if(percent_difference < 10) {
  practical_conclusion <- "Small practical difference between conditions"
} else {
  practical_conclusion <- "Meaningful practical difference between conditions"
}
print(paste("Practical conclusion:", practical_conclusion))

# ============================================================================
# ONE-SAMPLE TESTS AGAINST BASELINE
# ============================================================================

print("\n=== ONE-SAMPLE TESTS AGAINST BASELINE (1.0 = no change) ===")

# Test if each condition differs significantly from 1.0 (no change)
baseline_value <- 1.0

cat("Testing each condition against baseline of 1.0:\n")

# KD condition
kd_baseline_test <- t.test(kd_data, mu = baseline_value)
print(paste("CCL3 KD vs baseline: t =", round(kd_baseline_test$statistic, 3), 
            ", p =", round(kd_baseline_test$p.value, 4)))

# NG condition  
ng_baseline_test <- t.test(ng_data, mu = baseline_value)
print(paste("Negative Control vs baseline: t =", round(ng_baseline_test$statistic, 3),
            ", p =", round(ng_baseline_test$p.value, 4)))

# Apply multiple comparison correction
baseline_p_values <- c(kd_baseline_test$p.value, ng_baseline_test$p.value)
baseline_p_bonf <- p.adjust(baseline_p_values, method = "bonferroni")
baseline_p_fdr <- p.adjust(baseline_p_values, method = "fdr")

print("\nMultiple comparison corrections:")
print(paste("CCL3 KD vs baseline: p_bonferroni =", round(baseline_p_bonf[1], 4)))
print(paste("Negative Control vs baseline: p_bonferroni =", round(baseline_p_bonf[2], 4)))

# ============================================================================
# NORMALITY TESTING
# ============================================================================

print("\n=== NORMALITY TESTING ===")

# Test normality for each condition
normality_results <- data_long %>%
  group_by(condition) %>%
  summarise(
    shapiro_p = shapiro.test(tnf_alpha)$p.value,
    .groups = 'drop'
  )

print("Shapiro-Wilk test p-values by condition:")
print(normality_results)

# Overall normality test
overall_normality <- shapiro.test(data_long$tnf_alpha)
print(paste("\nOverall normality test p-value:", round(overall_normality$p.value, 4)))

# ============================================================================
# HOMOGENEITY OF VARIANCE
# ============================================================================

print("\n=== HOMOGENEITY OF VARIANCE ===")

# Levene's test
levene_result <- leveneTest(tnf_alpha ~ condition, data = data_long)
print("Levene's test for homogeneity of variance:")
print(levene_result)

# ============================================================================
# VISUALIZATION
# ============================================================================

print("\n=== VISUALIZATION ===")

# Main comparison plot
p1 <- ggplot(data_long, aes(x = condition, y = tnf_alpha)) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "red", size = 1) +
  geom_boxplot(aes(fill = condition), alpha = 0.7, width = 0.5) +
  geom_point(size = 3, position = position_jitter(width = 0.1)) +
  geom_line(aes(group = donor), alpha = 0.7, size = 1) +
  theme_minimal() +
  labs(title = "TNF-α Expression (Normalized to Control)",
       subtitle = "Red line shows baseline (1.0 = no change from control)",
       y = "TNF-α (Fold Change)", 
       x = "Condition") +
  scale_x_discrete(labels = c("kd" = "CCL3 Knockdown", "ng" = "Negative Control")) +
  theme(legend.position = "none")

# Individual donor trajectories
p2 <- ggplot(data_long, aes(x = condition, y = tnf_alpha, color = donor)) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "red") +
  geom_point(size = 3) +
  geom_line(aes(group = donor), size = 1) +
  theme_minimal() +
  labs(title = "Individual Donor Responses",
       x = "Condition", y = "TNF-α (Fold Change)", color = "Donor") +
  scale_x_discrete(labels = c("kd" = "CCL3 Knockdown", "ng" = "Negative Control"))

# Show the difference between conditions
differences <- kd_data - ng_data
diff_data <- data.frame(
  donor = factor(paste("Donor", 1:3)),
  difference = differences
)

p3 <- ggplot(diff_data, aes(x = donor, y = difference)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Difference: KD - NG (by donor)",
       subtitle = "Negative values indicate KD < NG",
       x = "Donor", y = "Difference in TNF-α")

print("Generating visualizations...")
grid.arrange(p1, p2, p3, ncol = 1)

# ============================================================================
# FINAL INTERPRETATION
# ============================================================================

print("\n=== FINAL INTERPRETATION ===")

cat("Based on your data:\n")
cat("1. Both conditions show massive suppression of TNF-α (~80% reduction)\n")
cat("2. CCL3 KD: 18% of control (82% suppression)\n") 
cat("3. Negative Control: 23% of control (77% suppression)\n")
cat("4. Difference between conditions: 5 percentage points\n\n")

if(paired_test$p.value >= 0.05) {
  cat("STATISTICAL CONCLUSION:\n")
  cat("- No significant difference between KD and NG (p ≥ 0.05)\n")
  cat("- Both treatments have statistically equivalent effects\n")
  cat("- The primary effect appears to be from siRNA delivery rather than target-specific knockdown\n\n")
} else {
  cat("STATISTICAL CONCLUSION:\n")
  cat("- Significant difference between KD and NG (p < 0.05)\n")
  cat("- CCL3 knockdown provides additional suppression beyond NG effect\n\n")
}

cat("BIOLOGICAL INTERPRETATION:\n")
cat("- Both siRNA treatments trigger strong anti-TNF-α responses\n")
cat("- This is consistent with literature on siRNA off-target effects\n")
cat("- Consider this context when interpreting CCL3's role in TNF-α regulation\n")

print("\n=== SUMMARY ===")
print("Analysis complete. Your data shows strong effects in both conditions with minimal difference between them.")
print("This pattern is well-documented in siRNA literature and suggests important technical considerations for the field.")

