# Select relevant columns
raw_il10_k <- select(raw, 1,3,6,7,10,12:22,24,25)

# 3. Read the sample metadata from tnfa_groups.csv
il10_groups_k <- read.csv('il10_groups_k.csv', header = TRUE, stringsAsFactors = FALSE)

# Set the row names to the gene IDs (first column), then remove that column.
rownames(raw_il10_k) <- raw_il10_k[[1]]
raw_il10_k <- raw_il10_k[, -1]

# Check that every column name in raw_tnfa (i.e. sample IDs) is found in tnfa_groups$id
if (all(colnames(raw_il10_k) %in% il10_groups_k$id)) {
  cat("All sample names in raw_ are present in _groups$id.\n")
} else {
  missing_cols <- setdiff(colnames(raw_il10_k), il10_groups_k$id)
  cat("The following sample names in raw_ are not found in _groups$id:\n")
  print(missing_cols)
}

# 4. Ensure that the 'group' column is a factor and set the reference level
il10_groups_k$group <- factor(il10_groups_k$group)
il10_groups_k$group <- relevel(il10_groups_k$group, ref = "low")

# 5. Create the DESeqDataSet from the raw counts.
#    Note: We do NOT use tidy=TRUE because our data is in wide format (genes x samples).
dds_10k <- DESeqDataSetFromMatrix(countData = raw_il10_k,
                              colData = il10_groups_k,
                              design = ~ group)

# Run the DESeq2 analysis (this will estimate size factors and perform normalization internally)
dds_10k <- DESeq(dds_10k)

# Retrieve the results (differential expression statistics)
res_10k <- results(dds_10k)
head(as.data.frame(res_10k))
summary(res_10k)

# Save the results to a CSV file
output_fn <- "C:/Users/brady/OneDrive/Desktop/TCM_RNA/Output/CorrelationTables/il10_pairwise_k.csv"
write.csv(as.data.frame(res), output_fn)

# Optional: Extract normalized counts (if you wish to inspect them later)
norm_counts <- counts(dds, normalized = TRUE)
head(norm_counts)

# Create a volcano plot using adjusted p-values
with(as.data.frame(res_10k), plot(log2FoldChange, -log10(padj),
                              pch = 20, main = "IL-10 (k means)"))

# Highlight significant points consistently using padj
with(subset(as.data.frame(res_10k), padj < 0.01 & log2FoldChange >= 2), 
     points(log2FoldChange, -log10(padj), pch = 20, col = "green"))
with(subset(as.data.frame(res_10k), padj < 0.01 & log2FoldChange <= -2), 
     points(log2FoldChange, -log10(padj), pch = 20, col = "red"))

#PCA plot
#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation

vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="group") #using the DESEQ2 plotPCA fxn we can