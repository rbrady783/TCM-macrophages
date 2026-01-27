
#Select relevant columns
raw_vegf_k <- select(raw, 1,3:8,14:17,19,20,22,25)

# 3. Read the sample metadata from vegf_k_groups.csv
vegf_k_groups <- read.csv('vegf_k_groups.csv', header = TRUE, stringsAsFactors = FALSE)

# Set the row names to the gene IDs (first column), then remove that column.
rownames(raw_vegf_k) <- raw_vegf_k[[1]]
raw_vegf_k <- raw_vegf_k[, -1]

# Check that every column name in raw_vegf_k (i.e. sample IDs) is found in vegf_k_groups$id
if (all(colnames(raw_vegf_k) %in% vegf_k_groups$id)) {
  cat("All sample names in raw_ are present in _groups$id.\n")
} else {
  missing_cols <- setdiff(colnames(raw_vegf_k), vegf_k_groups$id)
  cat("The following sample names in raw_ are not found in _groups$id:\n")
  print(missing_cols)
}

# 4. Ensure that the 'group' column is a factor and set the reference level
vegf_k_groups$group <- factor(vegf_k_groups$group)
vegf_k_groups$group <- relevel(vegf_k_groups$group, ref = "low")

# 5. Create the DESeqDataSet from the raw counts.
#    Note: We do NOT use tidy=TRUE because our data is in wide format (genes x samples).
dds_vegf_k <- DESeqDataSetFromMatrix(countData = raw_vegf_k,
                                   colData = vegf_k_groups,
                                   design = ~ group)

# Run the DESeq2 analysis (this will estimate size factors and perform normalization internally)
dds_vegf_k <- DESeq(dds_vegf_k)

# Retrieve the results (differential expression statistics)
res_vegf_k <- results(dds_vegf_k)
head(as.data.frame(res_vegf_k))
summary(res_vegf_k)

# Save the results to a CSV file
output_fn <- "C:/Users/brady/OneDrive/Desktop/TCM_RNA/Output/CorrelationTables/vegf_k_pairwise.csv"
write.csv(as.data.frame(res), output_fn)

# Optional: Extract normalized counts (if you wish to inspect them later)
norm_counts <- counts(dds, normalized = TRUE)
head(norm_counts)

# Create a volcano plot using adjusted p-values
with(as.data.frame(res_vegf_k), plot(log2FoldChange, -log10(padj),
                                   pch = 20, main = "VEGF (k means)"))

# Highlight significant points consistently using padj
with(subset(as.data.frame(res_vegf_k), padj < 0.01 & log2FoldChange >= 2), 
     points(log2FoldChange, -log10(padj), pch = 20, col = "green"))
with(subset(as.data.frame(res_vegf_k), padj < 0.01 & log2FoldChange <= -2), 
     points(log2FoldChange, -log10(padj), pch = 20, col = "red"))

#PCA plot
#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation

vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="group") #using the DESEQ2 plotPCA fxn we can