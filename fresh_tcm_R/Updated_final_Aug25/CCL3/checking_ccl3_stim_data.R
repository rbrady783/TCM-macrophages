# Load required packages
library(ggplot2)
library(dplyr)
library(car)        # For normality and homogeneity tests
#install.packages("moments")
library(moments)    # For skewness and kurtosis
library(gridExtra)  # For arranging plots
#install.packages("lmtest")
library(lmtest)

# Load your CSV data
data <- read.csv("recom_ccl3_tnfa.csv", stringsAsFactors = FALSE)
data$donor <- as.factor(data$donor)

cat("=============================================================\n")
cat("STEP 1: INITIAL DATA EXPLORATION\n")
cat("=============================================================\n")

# Basic data summary
str(data)
print(summary(data))

# Check for missing values
cat("\nMissing values:\n")
print(sapply(data, function(x) sum(is.na(x))))

# Check data ranges
cat("\nData ranges:\n")
cat("CCL3 dose range:", range(data$dose), "\n")
cat("TNF-α range:", range(data$tnfa), "\n")

cat("\n=============================================================\n")
cat("STEP 2: VISUALIZE RAW DATA DISTRIBUTION\n")
cat("=============================================================\n")

# Create multiple plots to assess data distribution
p1 <- ggplot(data, aes(x = tnfa)) +
  geom_histogram(bins = 8, fill = "lightblue", alpha = 0.7) +
  geom_density(aes(y = after_stat(count)), color = "red", linewidth = 1) +
  labs(title = "TNF-α Distribution", x = "TNF-α", y = "Count") +
  theme_minimal()

p2 <- ggplot(data, aes(x = dose)) +
  geom_histogram(bins = 8, fill = "lightgreen", alpha = 0.7) +
  labs(title = "CCL3 Dose Distribution", x = "CCL3 Dose", y = "Count") +
  theme_minimal()

p3 <- ggplot(data, aes(sample = tnfa)) +
  stat_qq() + stat_qq_line() +
  labs(title = "TNF-α Q-Q Plot") +
  theme_minimal()

p4 <- ggplot(data, aes(x = dose, y = tnfa)) +
  geom_point(aes(color = donor), size = 3) +
  geom_smooth(method = "loess", se = TRUE) +
  labs(title = "Raw Data: CCL3 vs TNF-α", 
       x = "CCL3 Dose", y = "TNF-α") +
  theme_minimal()

# Print plots
grid.arrange(p1, p2, p3, p4, ncol = 2)

cat("\n=============================================================\n")
cat("STEP 3: NORMALITY TESTS\n")
cat("=============================================================\n")

# Test normality of TNF-α
shapiro_tnfa <- shapiro.test(data$tnfa)
cat("TNF-α Shapiro-Wilk test:\n")
cat("  W =", round(shapiro_tnfa$statistic, 4), "\n")
cat("  p-value =", round(shapiro_tnfa$p.value, 4), "\n")
cat("  Interpretation:", ifelse(shapiro_tnfa$p.value < 0.05, 
                                "NOT normally distributed (p < 0.05)", 
                                "Appears normally distributed (p ≥ 0.05)"), "\n")

# Calculate skewness and kurtosis
tnfa_skew <- skewness(data$tnfa)
tnfa_kurt <- kurtosis(data$tnfa)

cat("\nTNF-α distribution characteristics:\n")
cat("  Skewness:", round(tnfa_skew, 3), 
    ifelse(abs(tnfa_skew) > 1, " (highly skewed)", 
           ifelse(abs(tnfa_skew) > 0.5, " (moderately skewed)", " (approximately symmetric)")), "\n")
cat("  Kurtosis:", round(tnfa_kurt, 3),
    ifelse(tnfa_kurt > 3, " (heavy-tailed)", 
           ifelse(tnfa_kurt < 3, " (light-tailed)", " (normal tails)")), "\n")

cat("\n=============================================================\n")
cat("STEP 4: TEST DIFFERENT TRANSFORMATIONS\n")
cat("=============================================================\n")

# Try common transformations
data$log_tnfa <- log(data$tnfa)
data$sqrt_tnfa <- sqrt(data$tnfa)
data$inv_tnfa <- 1/data$tnfa

# Test normality for each transformation
transformations <- list(
  "Original" = data$tnfa,
  "Log" = data$log_tnfa,
  "Square root" = data$sqrt_tnfa,
  "Inverse" = data$inv_tnfa
)

cat("Normality tests for different transformations:\n")
transformation_results <- data.frame(
  Transformation = character(),
  Shapiro_W = numeric(),
  Shapiro_p = numeric(),
  Skewness = numeric(),
  Normal = character(),
  stringsAsFactors = FALSE
)

for(i in 1:length(transformations)) {
  name <- names(transformations)[i]
  values <- transformations[[i]]
  
  # Skip if any infinite or NaN values
  if(any(!is.finite(values))) {
    cat(name, ": Contains non-finite values, skipping\n")
    next
  }
  
  shapiro_test <- shapiro.test(values)
  skew <- skewness(values)
  
  transformation_results <- rbind(transformation_results, 
                                  data.frame(
                                    Transformation = name,
                                    Shapiro_W = round(shapiro_test$statistic, 4),
                                    Shapiro_p = round(shapiro_test$p.value, 4),
                                    Skewness = round(skew, 3),
                                    Normal = ifelse(shapiro_test$p.value >= 0.05, "Yes", "No")
                                  ))
}

print(transformation_results)

# Find best transformation
best_transform <- transformation_results[which.max(transformation_results$Shapiro_p), ]
cat("\nBest transformation (highest Shapiro p-value):", best_transform$Transformation, "\n")

cat("\n=============================================================\n")
cat("STEP 5: VISUALIZE TRANSFORMATIONS\n")
cat("=============================================================\n")

# Create plots for all transformations
plot_list <- list()
for(i in 1:length(transformations)) {
  name <- names(transformations)[i]
  values <- transformations[[i]]
  
  if(any(!is.finite(values))) next
  
  plot_list[[name]] <- ggplot(data.frame(x = values), aes(x = x)) +
    geom_histogram(bins = 6, fill = "lightblue", alpha = 0.7) +
    geom_density(aes(y = after_stat(count)), color = "red") +
    labs(title = paste(name, "TNF-α"), x = paste(name, "TNF-α")) +
    theme_minimal()
}

if(length(plot_list) > 0) {
  do.call(grid.arrange, c(plot_list, ncol = 2))
}

cat("\n=============================================================\n")
cat("STEP 6: CHECK HOMOGENEITY OF VARIANCE\n")
cat("=============================================================\n")

# Test homogeneity of variance across dose groups
cat("Levene's test for homogeneity of variance:\n")

# Group doses for testing (since you have specific dose levels)
data$dose_group <- as.factor(data$dose)

levene_original <- leveneTest(tnfa ~ dose_group, data = data)
cat("Original TNF-α: F =", round(levene_original$`F value`[1], 4), 
    ", p =", round(levene_original$`Pr(>F)`[1], 4), "\n")

if(exists("log_tnfa", where = data)) {
  levene_log <- leveneTest(log_tnfa ~ dose_group, data = data)
  cat("Log TNF-α: F =", round(levene_log$`F value`[1], 4), 
      ", p =", round(levene_log$`Pr(>F)`[1], 4), "\n")
}

cat("\n=============================================================\n")
cat("STEP 7: CHECK LINEARITY ASSUMPTIONS\n")
cat("=============================================================\n")

# Plot residuals vs fitted for linearity check (simple linear model)
simple_lm <- lm(tnfa ~ dose, data = data)
data$residuals <- residuals(simple_lm)
data$fitted <- fitted(simple_lm)

p_linearity <- ggplot(data, aes(x = fitted, y = residuals)) +
  geom_point(size = 3) +
  geom_smooth(method = "loess", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Residuals vs Fitted (Check Linearity)",
       x = "Fitted Values", y = "Residuals") +
  theme_minimal()

print(p_linearity)

# Test for curvature
reset_test <- resettest(simple_lm, power = 2:3, type = "regressor")
cat("\nRESET test for linearity:\n")
cat("F =", round(reset_test$statistic, 4), ", p =", round(reset_test$p.value, 4), "\n")
cat("Interpretation:", ifelse(reset_test$p.value < 0.05, 
                              "Non-linear relationship detected", 
                              "Linear relationship appears adequate"), "\n")

cat("\n=============================================================\n")
cat("STEP 8: RECOMMENDATIONS\n")
cat("=============================================================\n")

cat("NORMALITY:\n")
if(best_transform$Normal == "Yes") {
  cat("✓ Use", best_transform$Transformation, "transformation for TNF-α\n")
} else {
  cat("⚠ No transformation achieves perfect normality\n")
  cat("Consider robust methods or nonparametric approaches\n")
}

cat("\nVARIANCE HOMOGENEITY:\n")
if(levene_original$`Pr(>F)`[1] >= 0.05) {
  cat("✓ Homogeneity assumption met for original data\n")
} else {
  cat("⚠ Variance heterogeneity detected - consider transformations or robust methods\n")
}

cat("\nLINEARITY:\n")
if(reset_test$p.value >= 0.05) {
  cat("✓ Linear relationship appears appropriate\n")
  cat("→ Proceed with linear mixed-effects models\n")
} else {
  cat("⚠ Non-linear relationship detected\n")
  cat("→ Consider polynomial or nonlinear mixed-effects models\n")
}

cat("\nFINAL RECOMMENDATION:\n")
if(best_transform$Transformation == "Original") {
  cat("Use original TNF-α values with standard mixed-effects models\n")
} else {
  cat("Use", best_transform$Transformation, "transformed TNF-α values\n")
}

cat("\nSAMPLE SIZE CONSIDERATION:\n")
cat("With n=16 observations (4 donors × 4 doses), you have limited power\n")
cat("Focus on effect sizes and confidence intervals, not just p-values\n")

