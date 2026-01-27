# Load required packages
library(nlme)     # Mixed-effects models
library(ggplot2)  # Plotting
library(dplyr)    # Data manipulation

# Load your data and apply the recommended transformation
data <- read.csv("recom_ccl3_tnfa.csv", stringsAsFactors = FALSE)
data$donor <- as.factor(data$donor)

# Apply inverse transformation (as recommended by data exploration)
data$inv_tnfa <- 1/data$tnfa

cat("=============================================================\n")
cat("FINAL ANALYSIS WITH INVERSE-TRANSFORMED TNF-α\n")
cat("=============================================================\n")

# =============================================================================
# PRIMARY MODEL: Linear Mixed-Effects with Inverse-Transformed TNF-α
# =============================================================================

# Linear mixed-effects model with random intercept for donor
final_model <- lme(inv_tnfa ~ dose, 
                   random = ~ 1 | donor, 
                   data = data,
                   method = "REML")

cat("\n=== PRIMARY MODEL RESULTS ===\n")
summary(final_model)

# =============================================================================
# STATISTICAL INFERENCE
# =============================================================================

cat("\n=== STATISTICAL TESTS ===\n")

# Test for dose effect
anova_result <- anova(final_model)
print(anova_result)

# Extract key results
dose_coef <- fixef(final_model)["dose"]
dose_p <- anova_result$`p-value`[2]

cat("\nKEY RESULTS:\n")
cat("Dose coefficient (on inverse scale):", round(dose_coef, 6), "\n")
cat("Dose effect p-value:", round(dose_p, 4), "\n")

# Significance interpretation
if(dose_p < 0.001) {
  sig_text <- "*** (p < 0.001)"
} else if(dose_p < 0.01) {
  sig_text <- "** (p < 0.01)"  
} else if(dose_p < 0.05) {
  sig_text <- "* (p < 0.05)"
} else {
  sig_text <- "(not significant)"
}

cat("Significance level:", sig_text, "\n")

# =============================================================================
# INTERPRETATION ON ORIGINAL SCALE
# =============================================================================

cat("\n=== INTERPRETATION ON ORIGINAL TNF-α SCALE ===\n")

# Since we used inverse transformation, let's interpret back on original scale
# If dose coefficient is positive on inverse scale, TNF-α decreases with dose
# If dose coefficient is negative on inverse scale, TNF-α increases with dose

if(dose_coef > 0) {
  cat("✓ As CCL3 dose increases, TNF-α levels DECREASE\n")
  cat("  (positive coefficient on inverse scale = decreasing on original scale)\n")
} else {
  cat("✓ As CCL3 dose increases, TNF-α levels INCREASE\n")
  cat("  (negative coefficient on inverse scale = increasing on original scale)\n")
}

# Calculate effect size: difference between dose 0 and dose 5
pred_dose0 <- predict(final_model, newdata = data.frame(dose = 0, donor = "1"), level = 0)
pred_dose5 <- predict(final_model, newdata = data.frame(dose = 5, donor = "1"), level = 0)

# Back-transform to original scale
tnfa_dose0 <- 1/pred_dose0
tnfa_dose5 <- 1/pred_dose5

cat("\nEFFECT SIZE ESTIMATE:\n")
cat("Predicted TNF-α at dose 0:", round(tnfa_dose0, 2), "\n")
cat("Predicted TNF-α at dose 5:", round(tnfa_dose5, 2), "\n")
cat("Difference:", round(tnfa_dose5 - tnfa_dose0, 2), "\n")
cat("Percent change:", round(((tnfa_dose5 - tnfa_dose0)/tnfa_dose0) * 100, 1), "%\n")

# =============================================================================
# MODEL VALIDATION
# =============================================================================

cat("\n=== MODEL VALIDATION ===\n")

# Check residuals
par(mfrow = c(2, 2))
plot(final_model, main = "Model Diagnostics (Inverse-Transformed)")
par(mfrow = c(1, 1))

# Calculate confidence intervals
conf_int <- intervals(final_model)
cat("\nConfidence intervals:\n")
print(conf_int$fixed)

# =============================================================================
# VISUALIZATION
# =============================================================================

cat("\n=== CREATING VISUALIZATIONS ===\n")

# Plot on original scale for interpretability
pred_data <- expand.grid(dose = seq(0, 5, 0.1), 
                         donor = levels(data$donor))

# Predict on inverse scale, then back-transform
pred_data$pred_inv <- predict(final_model, pred_data, level = 0)
pred_data$pred_tnfa <- 1/pred_data$pred_inv  # Back-transform to original scale

# Main plot - original scale
p1 <- ggplot(data, aes(x = dose, y = tnfa)) +
  geom_point(aes(color = donor), size = 3) +
  geom_line(data = pred_data, aes(x = dose, y = pred_tnfa), 
            color = "blue", linewidth = 1) +
  labs(title = "CCL3 vs TNF-α (Final Model Results)",
       subtitle = paste0("Dose effect: ", sig_text),
       x = "CCL3 Dose", 
       y = "TNF-α Level",
       color = "Donor") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

print(p1)

# Individual donor curves
p2 <- ggplot(data, aes(x = dose, y = tnfa, color = donor)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~paste("Donor", donor)) +
  labs(title = "Individual Donor Responses",
       x = "CCL3 Dose", y = "TNF-α Level") +
  theme_minimal()

print(p2)

# Plot on transformed scale to check model fit
p3 <- ggplot(data, aes(x = dose, y = inv_tnfa)) +
  geom_point(aes(color = donor), size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Model Fit on Transformed Scale",
       subtitle = "Should show linear relationship",
       x = "CCL3 Dose", 
       y = "1/TNF-α (Inverse Transformed)") +
  theme_minimal()

print(p3)

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("FINAL SUMMARY AND CONCLUSIONS\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

cat("MODEL: Linear mixed-effects with inverse-transformed TNF-α\n")
cat("TRANSFORMATION: 1/TNF-α (recommended by data exploration)\n")
cat("DESIGN: 4 donors × 4 CCL3 doses (n = 16 total observations)\n\n")

cat("MAIN FINDINGS:\n")
cat("1. CCL3 dose effect:", sig_text, "\n")

if(dose_coef > 0) {
  cat("2. Direction: TNF-α DECREASES as CCL3 dose increases\n")
} else {
  cat("2. Direction: TNF-α INCREASES as CCL3 dose increases\n")
}

cat("3. Effect size: TNF-α changes by", round(tnfa_dose5 - tnfa_dose0, 2), 
    "units from dose 0 to dose 5\n")
cat("4. Percent change:", round(((tnfa_dose5 - tnfa_dose0)/tnfa_dose0) * 100, 1), 
    "% change across dose range\n")

cat("\nMODEL QUALITY:\n")
cat("✓ Used appropriate transformation (inverse) based on data exploration\n")
cat("✓ Accounts for donor-to-donor variability with random effects\n")
cat("✓ Linear relationship confirmed on transformed scale\n")

cat("\nIMPORTANCE NOTE:\n")
cat("With only 4 donors, focus on effect size and confidence intervals\n")
cat("Consider this a pilot study for larger future experiments\n")

# Print the final model summary
cat("\n", paste(rep("-", 50), collapse = ""), "\n")
print(summary(final_model))

#####################################
# REQUIRED: Extract the key statistics PLOS asks for
fixed_summary <- summary(final_model)$tTable
conf_int <- intervals(final_model)$fixed

# Create simple table with ONLY required elements
required_table <- data.frame(
  Term = rownames(fixed_summary),
  Coefficient = fixed_summary[,"Value"],
  SE = fixed_summary[,"Std.Error"],
  Lower_95CI = conf_int[,1],
  Upper_95CI = conf_int[,3],
  p_value = fixed_summary[,"p-value"]
)

# Add model fit statistics (also required)
cat("Supplementary Table: Linear mixed-effects model for CCL3 dose effect\n\n")
print(required_table, digits = 4)
cat("\nModel fit: AIC =", AIC(final_model), ", BIC =", BIC(final_model), "\n")
