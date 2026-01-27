# Add Wald statistic, confidence intervals, and clean cytokine names
deg_data_complete <- deg_data %>%
  mutate(
    # Extract cytokine name from source_file
    cytokine = str_remove(source_file, "_quartile_results\\.csv"),
    
    # Clean up cytokine names for display
    cytokine = case_when(
      cytokine == "il.8" ~ "IL-8",
      cytokine == "il.10" ~ "IL-10",
      cytokine == "kc.like" ~ "KC-like",
      cytokine == "tnf.a" ~ "TNF-α",
      cytokine == "tgf.b" ~ "TGF-β",
      cytokine == "vegf" ~ "VEGF",
      cytokine == "ccl2" ~ "CCL2",
      TRUE ~ cytokine  # Keep as-is if not matched
    ),
    
    # Wald statistic = log2FC / SE
    stat = log2FoldChange / lfcSE,
    
    # 95% CI = estimate ± 1.96 × SE
    conf.low = log2FoldChange - (1.96 * lfcSE),
    conf.high = log2FoldChange + (1.96 * lfcSE)
  )

# Reorder columns for PLOS requirements
deg_data_complete <- deg_data_complete %>%
  select(cytokine, gene, Hugo_Symbol, baseMean, 
         log2FoldChange, lfcSE, conf.low, conf.high, 
         stat, pvalue, padj)

# Save complete table
write_csv(deg_data_complete, "S_Table_Complete_DEG_Results.csv")

# Summary
cat("Total DEGs:", nrow(deg_data_complete), "\n")
cat("Unique genes:", n_distinct(deg_data_complete$gene), "\n")
cat("\nDEGs per cytokine:\n")
print(table(deg_data_complete$cytokine))

library(tidyverse)

# Read the data
deg_data_complete <- read_csv("S_Table_Complete_DEG_Results.csv")

# Remove rows where all key statistical columns are NA
deg_data_clean <- deg_data_complete %>%
  filter(!(is.na(log2FoldChange) & is.na(pvalue) & is.na(padj)))

# Or more comprehensive - remove if ANY of the key columns are NA
deg_data_clean <- deg_data_complete %>%
  filter(!is.na(log2FoldChange) & !is.na(pvalue) & !is.na(padj))

# Check how many rows were removed
cat("Original rows:", nrow(deg_data_complete), "\n")
cat("After removing NAs:", nrow(deg_data_clean), "\n")
cat("Rows removed:", nrow(deg_data_complete) - nrow(deg_data_clean), "\n")

# Save clean version
write_csv(deg_data_clean, "S_Table_Complete_DEG_Results_Clean.csv")
