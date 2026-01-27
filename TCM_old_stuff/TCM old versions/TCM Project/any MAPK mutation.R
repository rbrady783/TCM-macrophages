# Define the columns of interest
cols <- c("BRAF", "KRAS", "NRAS", "ERBB2", "PTPN11", "NF1", "KIT")

# Create the new column by checking each row:
joined_data$any_mutation <- apply(joined_data[, cols], 1, function(x) {
  # Standardize the values to lower case and trim any extra whitespace:
  x_clean <- tolower(trimws(x))
  # If any of the values is "yes", return 1; otherwise return 0.
  if (any(x_clean == "yes", na.rm = TRUE)) {
    return(1)
  } else {
    return(0)
  }
})

# Convert the any_mutation column to a factor with labels "No" for 0 and "Yes" for 1
joined_data$any_mutation <- factor(joined_data$any_mutation, levels = c(0, 1),
                                   labels = c("No", "Yes"))

# Get the names of all numeric columns
numeric_vars <- names(joined_data)[sapply(joined_data, is.numeric)]
# Optionally, inspect the numeric variables
print(numeric_vars)

library(ggplot2)

for (num in numeric_vars) {
  # Create a boxplot with jittered points comparing the numeric variable by any_mutation
  p <- ggplot(joined_data, aes(x = any_mutation, y = .data[[num]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
    labs(title = paste("Comparison of", num, "by Any Mutation"),
         x = "Any Mutation", y = num) +
    theme_minimal()
  
  print(p)
  
  # Perform a t-test comparing the numeric variable between the two groups
  test_result <- t.test(joined_data[[num]] ~ joined_data$any_mutation)
  cat("T-test comparing", num, "by Any Mutation:\n")
  print(test_result)
  cat("--------------------------------------------------\n")
}
