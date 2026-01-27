library(tidyverse)
library(readxl)

# Get all sheet names (7 cytokines)
sheet_names <- excel_sheets("Combined_GSEA_Results_2025-08-29.xlsx")
print(sheet_names)

# Read all sheets and combine
gsea_all <- map_df(sheet_names, ~ {
  read_excel("Combined_GSEA_Results_2025-08-29.xlsx", sheet = .x)
})

# Filter for enriched gene sets (qvalue ≤ 0.25 and |NES| ≥ 1.0)
enriched <- gsea_all %>%
  filter(qvalue <= 0.25 & abs(NES) >= 1.0)

# Total enriched gene sets
total_enriched <- nrow(enriched)
cat("Total enriched gene sets:", total_enriched, "\n\n")

# Count by cytokine
enriched_by_cytokine <- enriched %>%
  group_by(Cytokine) %>%
  summarise(n_gene_sets = n()) %>%
  arrange(desc(n_gene_sets))

print(enriched_by_cytokine)

# Range per cytokine
min_per_cytokine <- min(enriched_by_cytokine$n_gene_sets)
max_per_cytokine <- max(enriched_by_cytokine$n_gene_sets)
cat("\nRange per cytokine:", min_per_cytokine, "-", max_per_cytokine, "\n")

# Breakdown by collection
enriched_by_collection <- enriched %>%
  group_by(Collection) %>%
  summarise(n_gene_sets = n())

cat("\nEnriched gene sets by collection:\n")
print(enriched_by_collection)

