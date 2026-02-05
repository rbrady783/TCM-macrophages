# Identify numeric and factor variables in joined_data
num_vars <- names(joined_data)[sapply(joined_data, is.numeric)]
fac_vars <- names(joined_data)[sapply(joined_data, is.factor)]

# Open a sink to capture output to a text file
sink("analysis_results_filtered.txt")

for(n in num_vars) {
  for(f in fac_vars) {
    
    cat("-----------------------------------------------------\n")
    cat("Analysis for", n, "by", f, "\n")
    
    # Print group counts for the factor variable
    cat("Group counts for", f, ":\n")
    print(table(joined_data[[f]]))
    
    # Get the counts of observations per factor level
    group_counts <- table(joined_data[[f]])
    
    # Only run tests if every group has at least 2 observations
    if(all(group_counts >= 2)) {
      if(length(levels(joined_data[[f]])) == 2) {
        # Perform a t-test for two-group comparisons
        test_result <- t.test(joined_data[[n]] ~ joined_data[[f]])
        p_val <- test_result$p.value
        if (p_val < 0.1) {
          cat("\nT-test comparing", n, "by", f, "(p-value:", p_val, "):\n")
          print(test_result)
        }
      } else if(length(levels(joined_data[[f]])) > 2) {
        # Perform an ANOVA for factors with more than two levels
        model <- aov(joined_data[[n]] ~ joined_data[[f]])
        anova_summary <- summary(model)[[1]]
        # Extract the p-value from the first row of the summary table
        p_val <- anova_summary$`Pr(>F)`[1]
        if (p_val < 0.1) {
          cat("\nANOVA comparing", n, "by", f, "(p-value:", p_val, "):\n")
          print(anova_summary)
        }
      }
    } else {
      cat("\nSkipping test for", n, "by", f, ": not enough observations in one or more groups.\n")
    }
    
    cat("-----------------------------------------------------\n\n")
  }
}

# Close the sink so that output returns to the console
sink()
