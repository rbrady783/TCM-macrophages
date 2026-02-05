
#Tumor type versus secretory

joined_data <- joined_data %>% rename(tumor.type = Tumor.type)

# Ensure tumor.type is a factor
joined_data$tumor.type <- as.factor(joined_data$tumor.type)

# Identify numeric columns in joined_data
numeric_vars <- names(joined_data)[sapply(joined_data, is.numeric)]

# =============================
# A. Save Textual Output to a File
# =============================

# Open a sink to capture all printed output to a text file
sink("tumor_type_analysis.txt")

for(num in numeric_vars) {
  
  cat("-----------------------------------------------------\n")
  cat("Analysis for", num, "by tumor.type:\n")
  
  # Print group counts for tumor.type
  cat("Group counts for tumor.type:\n")
  print(table(joined_data$tumor.type))
  
  # Run the statistical test
  if(length(levels(joined_data$tumor.type)) == 2) {
    test_result <- t.test(joined_data[[num]] ~ joined_data$tumor.type)
    cat("\nT-test comparing", num, "by tumor.type:\n")
    print(test_result)
  } else {
    model <- aov(joined_data[[num]] ~ joined_data$tumor.type)
    cat("\nANOVA comparing", num, "by tumor.type:\n")
    print(summary(model))
  }
  
  cat("-----------------------------------------------------\n\n")
}

# Close the sink so that output returns to the console
sink()

# =============================
# B. Save Plots to a PDF
# =============================

# Open a PDF device to save plots
pdf("tumor_type_plots.pdf", width = 8, height = 6)

for(num in numeric_vars) {
  p <- ggplot(joined_data, aes(x = tumor.type, y = .data[[num]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
    labs(title = paste("Comparison of", num, "across Tumor Types"),
         x = "Tumor Type", y = num) +
    theme_minimal()
  print(p)
}

dev.off()
