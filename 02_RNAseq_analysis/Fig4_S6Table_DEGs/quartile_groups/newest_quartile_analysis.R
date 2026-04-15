# ==========================================
# QUARTILE ANALYSIS FOR CYTOKINES
# ==========================================
# Produces two output tables:
#   1. n=25 cell lines (all available) -> S6 Table in manuscript
#   2. n=23 cell lines (Vogel and CMT27 excluded -- no RNA-seq data)
#      -> input for downstream DEG analysis in 02_RNAseq_analysis
# ==========================================

library(dplyr)
library(tidyr)

# Read data (contains 25 cell lines)
dat <- read.csv("modz_no_ctrls_histo_mean.csv", stringsAsFactors = FALSE)

# ------------------------------------------
# HELPER FUNCTIONS (unchanged from original)
# ------------------------------------------

# Determine which statistical test to use based on normality and variance
select_test <- function(top_values, bottom_values) {
  if (length(top_values) >= 3 && length(bottom_values) >= 3) {
    shapiro_top <- shapiro.test(top_values)
    shapiro_bottom <- shapiro.test(bottom_values)
    
    if (shapiro_top$p.value > 0.05 && shapiro_bottom$p.value > 0.05) {
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
    return("Mann-Whitney U test")
  }
}

# Run the appropriate test
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

# ------------------------------------------
# MAIN ANALYSIS FUNCTION
# ------------------------------------------
# Uses RANK-BASED selection (top 6 / bottom 6 cell lines per cytokine)
# to match the published S6 Table methodology, rather than a true 
# quantile-based cutoff.

run_quartile_analysis <- function(input_data, n_per_group = 6) {
  results_list <- list()
  
  for (cyt in unique(input_data$cytokine)) {
    cyt_data <- input_data %>% 
      filter(cytokine == cyt) %>%
      arrange(mean)
    
    # Rank-based selection: bottom 6 and top 6 cell lines
    bottom_quartile <- head(cyt_data, n_per_group)
    top_quartile <- tail(cyt_data, n_per_group)
    
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
  
  bind_rows(results_list)
}

# ------------------------------------------
# RUN 1: ALL 25 CELL LINES (S6 Table)
# ------------------------------------------
cat("=== Running analysis on all 25 cell lines (for S6 Table) ===\n")
results_n25 <- run_quartile_analysis(dat)
print(results_n25)
write.csv(results_n25, "S6_Table_quartile_groups_n25.csv", row.names = FALSE)

# ------------------------------------------
# RUN 2: 23 CELL LINES (DEG analysis subset)
# Vogel and CMT27 had no RNA-seq data, so were excluded from DEG analysis
# ------------------------------------------
cat("\n=== Running analysis on 23 cell lines (excluding Vogel and CMT27) ===\n")
dat_n23 <- dat %>%
  filter(!grepl("Vogel|CMT27", treatment, ignore.case = TRUE))

results_n23 <- run_quartile_analysis(dat_n23)
print(results_n23)
write.csv(results_n23, "quartile_groups_n23_for_DEG_analysis.csv", row.names = FALSE)

# ------------------------------------------
# SUMMARY
# ------------------------------------------
cat("\n=== Analysis complete ===\n")
cat(sprintf("n=25: %d cytokines, %d significant (p < 0.05)\n",
            nrow(results_n25),
            sum(results_n25$p_value < 0.05, na.rm = TRUE)))
cat(sprintf("n=23: %d cytokines, %d significant (p < 0.05)\n",
            nrow(results_n23),
            sum(results_n23$p_value < 0.05, na.rm = TRUE)))

