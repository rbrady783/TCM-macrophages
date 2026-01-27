library(tidyverse)

# Read your DEG file
deg_data <- read_csv("all_DEGS.csv")

# Filter for significant DEGs
# Adjusted p ≤ 0.05 (you said "greater than" but I think you meant "less than" for significance)
# |log2FC| ≥ 1 (i.e., log2FC ≥ 1 OR log2FC ≤ -1)
significant_degs <- deg_data %>%
  filter(padj <= 0.05 & (log2FoldChange >= 1 | log2FoldChange <= -1))

# Total number of DEGs (may include same gene in multiple cytokines)
total_degs <- nrow(significant_degs)
cat("Total significant DEGs:", total_degs, "\n")

# Number of UNIQUE genes (if you have gene IDs)
# Adjust column name to match your data (e.g., "gene_id", "ensembl_id", "GeneID")
unique_genes <- n_distinct(significant_degs$gene)  # Change "gene_id" to your column name
cat("Unique genes:", unique_genes, "\n")

# Count by cytokine (if you have a cytokine column)
degs_by_cytokine <- significant_degs %>%
  group_by(source_file) %>%  # Change "cytokine" to your column name if different
  summarise(n_DEGs = n()) %>%
  arrange(desc(n_DEGs))

print(degs_by_cytokine)

# Count up vs down-regulated
up_down_summary <- significant_degs %>%
  mutate(direction = ifelse(log2FoldChange > 0, "Up", "Down")) %>%
  count(direction)

print(up_down_summary)

###################################

library(tidyverse)

# Assuming your DEG data is already filtered for significance
# (padj <= 0.05 and |log2FC| >= 1)

# Count how many cytokines each gene is significant for
gene_frequency <- significant_degs %>%
  group_by(gene) %>%  # Change "gene_id" to your gene column name
  summarise(
    n_cytokines = n_distinct(source_file),  # Change "cytokine" to your column name
    cytokines = paste(unique(source_file), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cytokines))

# Show genes expressed in all 7 cytokines
genes_in_all_7 <- gene_frequency %>%
  filter(n_cytokines == 7)

cat("Genes significant for all 7 cytokines:", nrow(genes_in_all_7), "\n\n")
print(genes_in_all_7)

# Show distribution (how many genes in 1, 2, 3... cytokines)
frequency_dist <- gene_frequency %>%
  count(n_cytokines) %>%
  rename(genes_count = n)

cat("\nDistribution of genes across cytokines:\n")
print(frequency_dist)

# Genes in 6 or more cytokines
genes_in_many <- gene_frequency %>%
  filter(n_cytokines >= 6)

cat("\nGenes significant for 6+ cytokines:\n")
print(genes_in_many)

