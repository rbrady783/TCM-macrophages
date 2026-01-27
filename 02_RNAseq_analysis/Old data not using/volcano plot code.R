library(ggplot2)

# Extract normalized counts from your DESeqDataSet (dds_tnf)
norm_counts <- counts(dds_tnf, normalized = TRUE)
norm_df <- as.data.frame(norm_counts)
norm_df$gene_id <- rownames(norm_df)

# Specify the gene of interest by its ENSCAF ID
gene <- "ENSCAFG00000029471"  # Replace with your specific ENSCAF ID

# Subset the data for the gene of interest and reshape from wide to long format
gene_data <- norm_df %>%
  filter(gene_id == gene) %>%
  pivot_longer(cols = -gene_id, names_to = "sample", values_to = "normalized_count")

# Merge with sample metadata (tnfa_groups); assumes tnfa_groups has a column 'id' matching the sample names
gene_data <- merge(gene_data, tnfa_groups, by.x = "sample", by.y = "id")

# Create a box-and-whisker plot with jittered individual points
ggplot(gene_data, aes(x = group, y = normalized_count)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue") +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  labs(title = paste("Normalized Expression for CLIC1"),
       x = "Group (high vs low)",
       y = "Normalized Counts") +
  theme_classic()

