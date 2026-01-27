# Quartile Comparison Table with Proper Statistical Tests
# Purpose: Create publication table showing cell line stimulators and appropriate statistical tests

# Load required libraries
library(dplyr)
library(readr)
library(knitr)
library(kableExtra)

# Function to get top and bottom quartile stimulators (6 samples each)
get_quartile_stimulators <- function(means_file, quartile_size = 6) {
  
  # Read the means data
  means_data <- read_csv(means_file, show_col_types = FALSE)
  
  cat("Processing cytokine stimulators from treatment data...\n")
  
  # Get cytokine columns (exclude treatment column)
  cytokine_cols <- colnames(means_data)[colnames(means_data) != "treatment"]
  
  results_list <- list()
  
  for(cytokine in cytokine_cols) {
    
    cat(paste("Analyzing", cytokine, "quartiles...\n"))
    
    # Sort treatments by cytokine value
    cytokine_data <- means_data %>%
      select(treatment, !!sym(cytokine)) %>%
      arrange(desc(!!sym(cytokine)))
    
    # Get top quartile (highest 6 values)
    top_quartile <- cytokine_data %>%
      head(quartile_size) %>%
      pull(treatment)
    
    # Get bottom quartile (lowest 6 values)
    bottom_quartile <- cytokine_data %>%
      tail(quartile_size) %>%
      arrange(!!sym(cytokine)) %>%  # Sort ascending for bottom
      pull(treatment)
    
    results_list[[cytokine]] <- list(
      cytokine = cytokine,
      top_stimulators = paste(top_quartile, collapse = ", "),
      bottom_stimulators = paste(bottom_quartile, collapse = ", "),
      top_values = cytokine_data %>% head(quartile_size) %>% pull(!!sym(cytokine)),
      bottom_values = cytokine_data %>% tail(quartile_size) %>% pull(!!sym(cytokine))
    )
    
    cat(paste("  - Top quartile (n=6):", paste(top_quartile, collapse = ", "), "\n"))
    cat(paste("  - Bottom quartile (n=6):", paste(bottom_quartile, collapse = ", "), "\n"))
  }
  
  return(results_list)
}

# Function to determine appropriate statistical test and calculate p-value
perform_statistical_test <- function(top_values, bottom_values) {
  
  # Check sample sizes
  n_top <- length(top_values)
  n_bottom <- length(bottom_values)
  
  if(n_top < 3 || n_bottom < 3) {
    return(list(test = "Insufficient data", p_value = 1.0))
  }
  
  # Test for normality using Shapiro-Wilk (if n >= 3 and n <= 50)
  if(n_top >= 3 && n_top <= 50 && n_bottom >= 3 && n_bottom <= 50) {
    shapiro_top <- shapiro.test(top_values)
    shapiro_bottom <- shapiro.test(bottom_values)
    
    top_normal <- shapiro_top$p.value > 0.05
    bottom_normal <- shapiro_bottom$p.value > 0.05
  } else {
    # For larger samples, could use other normality tests
    top_normal <- TRUE
    bottom_normal <- TRUE
  }
  
  # Test for equal variances (F-test)
  var_test <- var.test(top_values, bottom_values)
  equal_variances <- var_test$p.value > 0.05
  
  # Choose appropriate test
  if(top_normal && bottom_normal) {
    if(equal_variances) {
      # Student's t-test (equal variances)
      test_result <- t.test(top_values, bottom_values, var.equal = TRUE)
      test_name <- "Student's t-test"
    } else {
      # Welch's t-test (unequal variances) 
      test_result <- t.test(top_values, bottom_values, var.equal = FALSE)
      test_name <- "Welch's t-test"
    }
  } else {
    # Mann-Whitney U test (non-parametric)
    test_result <- wilcox.test(top_values, bottom_values)
    test_name <- "Mann-Whitney"
  }
  
  return(list(test = test_name, p_value = test_result$p.value))
}

# Function to create complete quartile comparison table
create_quartile_comparison_table <- function(means_file = "means_modz.csv", quartile_size = 6) {
  
  # Get stimulators from the means data
  stimulator_data <- get_quartile_stimulators(means_file, quartile_size)
  
  results_list <- list()
  
  for(cytokine_col in names(stimulator_data)) {
    
    stim_info <- stimulator_data[[cytokine_col]]
    
    # Convert cytokine column name to display name
    display_name <- case_when(
      cytokine_col == "TNF-a" ~ "TNF-α",
      cytokine_col == "TGF-b" ~ "TGF-β", 
      TRUE ~ cytokine_col
    )
    
    cat(paste("\nPerforming statistical test for", display_name, "...\n"))
    
    # Perform statistical test between top and bottom quartiles
    stat_result <- perform_statistical_test(stim_info$top_values, stim_info$bottom_values)
    
    cat(paste("  - Test used:", stat_result$test, "\n"))
    cat(paste("  - p-value:", format(stat_result$p_value, scientific = TRUE), "\n"))
    
    results_list[[display_name]] <- data.frame(
      Cytokine = display_name,
      Top_Stimulators = stim_info$top_stimulators,
      Bottom_Stimulators = stim_info$bottom_stimulators,
      Test = stat_result$test,
      p_value = stat_result$p_value,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine all results
  final_table <- do.call(rbind, results_list)
  
  # Format p-values for display
  final_table$p_value_formatted <- sapply(final_table$p_value, function(p) {
    if(is.na(p) || p >= 1) {
      "NS"
    } else if(p < 0.001) {
      "<0.001"
    } else if(p < 0.01) {
      sprintf("%.4f", p)
    } else {
      sprintf("%.3f", p)
    }
  })
  
  return(final_table)
}

# Function to create publication-ready table
create_publication_table <- function(summary_df, save_files = TRUE) {
  
  # Select and rename columns for final table
  table_final <- summary_df %>%
    select(Top_Stimulators, Bottom_Stimulators, Test, p_value_formatted) %>%
    mutate(Cytokine = summary_df$Cytokine) %>%
    select(Cytokine, Top_Stimulators, Bottom_Stimulators, Test, p_value_formatted) %>%
    rename(
      "Cytokine" = "Cytokine",
      "Top Stimulators" = "Top_Stimulators",
      "Bottom Stimulators" = "Bottom_Stimulators", 
      "Test" = "Test",
      "p-value" = "p_value_formatted"
    )
  
  # Print simple table to console
  cat("\n=== QUARTILE COMPARISON SUMMARY ===\n")
  print(kable(table_final, format = "simple"))
  
  if(save_files) {
    # Save as CSV for easy editing
    write_csv(table_final, "quartile_comparison_table.csv")
    
    # Create HTML table with borders but no color formatting
    simple_html <- kable(table_final, 
                         format = "html",
                         caption = "Statistical Comparison of Top vs Bottom Quartile Cytokine-Stimulating Treatments",
                         table.attr = 'border="1" cellpadding="5" cellspacing="0" style="border-collapse: collapse;"')
    
    # Save HTML table
    writeLines(as.character(simple_html), "quartile_comparison_table.html")
    
    cat("\nFiles saved:\n")
    cat("  - quartile_comparison_table.csv\n")
    cat("  - quartile_comparison_table.html\n")
  }
  
  return(table_final)
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

# Generate the complete quartile comparison table
cat("Creating quartile comparison table with proper statistical tests...\n")
summary_data <- create_quartile_comparison_table("means_modz.csv", quartile_size = 6)

# Create the publication table
final_table <- create_publication_table(summary_data)

# Print methodology summary
cat("\n=== METHODOLOGY SUMMARY ===\n")
cat("Quartile Definition: Top 6 and bottom 6 treatments based on cytokine expression\n")
cat("Statistical Tests Used:\n")
cat("  - Student's t-test: Normal data, equal variances\n")
cat("  - Welch's t-test: Normal data, unequal variances\n") 
cat("  - Mann-Whitney: Non-normal data\n")
cat("Test Selection: Based on Shapiro-Wilk normality test and F-test for equal variances\n")
