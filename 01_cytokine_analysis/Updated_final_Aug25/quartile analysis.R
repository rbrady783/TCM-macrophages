# ==========================================
# QUARTILE ANALYSIS FOR CYTOKINES
# ==========================================

library(dplyr)
library(tidyr)

# Read data
dat <- read.csv("modz_no_ctrls_histo_mean.csv", stringsAsFactors = FALSE)

# Function to determine which statistical test to use
select_test <- function(top_values, bottom_values) {
  # Check normality with Shapiro-Wilk test
  if (length(top_values) >= 3 && length(bottom_values) >= 3) {
    shapiro_top <- shapiro.test(top_values)
    shapiro_bottom <- shapiro.test(bottom_values)
    
    # If both are normal (p > 0.05), check equal variances
    if (shapiro_top$p.value > 0.05 && shapiro_bottom$p.value > 0.05) {
      # Check equal variances with F-test
      var_test <- var.test(top_values, bottom_values)
      
      if (var_test$p.value > 0.05) {
        return("Student's t-test")
      } else {
        return("Welch's t-test")
      }
    } else {
      return("Mann-Whitney U test")
    }
  } else {
    # Too few samples, use non-parametric
    return("Mann-Whitney U test")
  }
}

# Function to run the appropriate test
run_test <- function(top_values, bottom_values, test_name) {
  if (test_name == "Student's t-test") {
    result <- t.test(top_values, bottom_values, var.equal = TRUE)
  } else if (test_name == "Welch's t-test") {
    result <- t.test(top_values, bottom_values, var.equal = FALSE)
  } else {
    result <- wilcox.test(top_values, bottom_values)
  }
  return(result$p.value)
}

# Analyze each cytokine
results_list <- list()

for (cyt in unique(dat$cytokine)) {
  # Filter for this cytokine and remove Vogel/CMT27
  cyt_data <- dat %>%
    filter(cytokine == cyt) %>%
    filter(!grepl("Vogel|CMT27", treatment, ignore.case = TRUE))
  
  # Calculate quartiles
  q1 <- quantile(cyt_data$mean, 0.25, na.rm = TRUE)
  q3 <- quantile(cyt_data$mean, 0.75, na.rm = TRUE)
  
  # Get top and bottom quartile data
  bottom_quartile <- cyt_data %>% filter(mean <= q1)
  top_quartile <- cyt_data %>% filter(mean >= q3)
  
  # Get cell lines (treatments)
  bottom_lines <- paste(bottom_quartile$treatment, collapse = "; ")
  top_lines <- paste(top_quartile$treatment, collapse = "; ")
  
  # Get values for testing
  bottom_values <- bottom_quartile$mean
  top_values <- top_quartile$mean
  
  # Determine and run test
  if (length(top_values) >= 2 && length(bottom_values) >= 2) {
    test_used <- select_test(top_values, bottom_values)
    p_value <- run_test(top_values, bottom_values, test_used)
  } else {
    test_used <- "Not enough data"
    p_value <- NA
  }
  
  # Store results
  results_list[[length(results_list) + 1]] <- data.frame(
    cytokine = cyt,
    bottom_quartile_cell_lines = bottom_lines,
    top_quartile_cell_lines = top_lines,
    n_bottom = length(bottom_values),
    n_top = length(top_values),
    test_used = test_used,
    p_value = p_value,
    stringsAsFactors = FALSE
  )
}

# Combine results
results_df <- bind_rows(results_list)

# Print results
print(results_df)

# Save to CSV
write.csv(results_df, "quartile_analysis_results.csv", row.names = FALSE)

cat("\nAnalysis complete! Results saved to 'quartile_analysis_results.csv'\n")
cat("\nSummary:\n")
cat(sprintf("  Total cytokines analyzed: %d\n", nrow(results_df)))
cat(sprintf("  Significant differences (p < 0.05): %d\n", 
            sum(results_df$p_value < 0.05, na.rm = TRUE)))
