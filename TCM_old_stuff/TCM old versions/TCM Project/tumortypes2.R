# 1. Convert the new column to a factor with desired labels.
joined_data$histo <- factor(joined_data$histo, 
                                  levels = c(1, 2, 3, 4),
                                  labels = c("sarcoma", "melanoma", "hematop", "carcinoma"))

# 2. Identify numeric columns in the joined_data data frame.
num_vars <- names(joined_data)[sapply(joined_data, is.numeric)]

# Optionally, check which numeric variables were found:
cat("Numeric columns:\n")
print(num_vars)

# 3. Compare each numeric variable with cancer_type.
library(ggplot2)

for(num in num_vars) {
  # Create a boxplot with jittered points.
  p <- ggplot(joined_data, aes(x = histo, y = .data[[num]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
    labs(title = paste("Comparison of", num, "across Cancer Types"),
         x = "Cancer Type", y = num) +
    theme_minimal()
  
  # Print the plot.
  print(p)
  
  # Since cancer_type has more than 2 levels, run an ANOVA.
  model <- aov(joined_data[[num]] ~ joined_data$histo)
  cat("\nANOVA comparing", num, "by histo:\n")
  print(summary(model))
  cat("-----------------------------------------------------\n")
}

#Print plots
# Open a PDF device to save the plots
pdf("histo.pdf", width = 8, height = 6)

# Loop over each numeric variable
for(num in num_vars) {
  p <- ggplot(joined_data, aes(x = histo, y = .data[[num]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
    labs(title = paste("Comparison of", num, "across Cancer Types"),
         x = "Cancer Type", y = num) +
    theme_minimal()
  
  print(p)  # This sends the plot to the PDF device
}

# Close the PDF device so the file is finalized
dev.off()
