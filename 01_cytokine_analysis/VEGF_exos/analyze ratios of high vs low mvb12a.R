# Your data (all 8 lines)
high_ratio <- c(17.723250130, 7.013631326, 14.357836050, 16.949359520)  
low_ratio <- c(5.775468997, 3.215865375, 5.840611760, 9.957062762)

# OPTION 1: Direct ratio (untransformed)
cat("=== OPTION 1: Untransformed Ratios ===\n")
high_mean <- mean(high_ratio)
high_sem <- sd(high_ratio) / sqrt(length(high_ratio))
low_mean <- mean(low_ratio)
low_sem <- sd(low_ratio) / sqrt(length(low_ratio))

cat("High: ", round(high_mean, 2), "±", round(high_sem, 2), "\n")
cat("Low: ", round(low_mean, 2), "±", round(low_sem, 2), "\n")

# Check normality
shapiro.test(high_ratio)
shapiro.test(low_ratio)

# Test
test1 <- t.test(high_ratio, low_ratio, var.equal = TRUE)
print(test1)

# OPTION 2: Log-transformed ratios
cat("\n=== OPTION 2: Log-Transformed Ratios ===\n")
high_log <- log(high_ratio)
low_log <- log(low_ratio)

high_log_mean <- mean(high_log)
high_log_sem <- sd(high_log) / sqrt(length(high_log))
low_log_mean <- mean(low_log)
low_log_sem <- sd(low_log) / sqrt(length(low_log))

cat("High (log): ", round(high_log_mean, 2), "±", round(high_log_sem, 2), "\n")
cat("Low (log): ", round(low_log_mean, 2), "±", round(low_log_sem, 2), "\n")

# Check normality on log scale
shapiro.test(high_log)
shapiro.test(low_log)

# Test on log scale
test2 <- t.test(high_log, low_log, var.equal = TRUE)
print(test2)

# Back-transform for interpretation
cat("\nGeometric mean ratio High:", round(exp(high_log_mean), 2), "\n")
cat("Geometric mean ratio Low:", round(exp(low_log_mean), 2), "\n")

