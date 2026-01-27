
#Mutations versus secretory products

# Identify numeric and factor variables in joined_data
num_vars <- names(joined_data)[sapply(joined_data, is.numeric)]
fac_vars  <- names(joined_data)[sapply(joined_data, is.factor)]

# Open a sink to capture all printed output to a text file
sink("analysis_results2.txt")

# Loop over each numeric variable
for(n in num_vars) {
  # Loop over each factor variable
  for(f in fac_vars) {
    
    cat("-----------------------------------------------------\n")
    cat("Analysis for", n, "by", f, "\n")
    
    # Print group counts for the current factor variable
    cat("Group counts for", f, ":\n")
    print(table(joined_data[[f]]))
    
    # Get the counts of observations per factor level
    group_counts <- table(joined_data[[f]])
    
    # Only run tests if every group has at least 2 observations
    if(all(group_counts >= 2)) {
      if(length(levels(joined_data[[f]])) == 2) {
        test_result <- t.test(joined_data[[n]] ~ joined_data[[f]])
        cat("\nT-test comparing", n, "by", f, ":\n")
        print(test_result)
      } else if(length(levels(joined_data[[f]])) > 2) {
        model <- aov(joined_data[[n]] ~ joined_data[[f]])
        cat("\nANOVA comparing", n, "by", f, ":\n")
        print(summary(model))
      }
    } else {
      cat("\nSkipping test for", n, "by", f, 
          ": not enough observations in one or more groups.\n")
    }
    
    cat("-----------------------------------------------------\n\n")
  }
}

# Close the sink so that output returns to the console
sink()
