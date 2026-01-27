# Select relevant columns
raw_il8_k <- select(raw, 1,3:5,9,16:23,25)

# 3. Read the sample metadata from il8_k_groups.csv
il8_k_groups <- read.csv('il8_k_groups.csv', header = TRUE, stringsAsFactors = FALSE)

# Set the row names to the gene IDs (first column), then remove that column.
rownames(raw_il8_k) <- raw_il8_k[[1]]
raw_il8_k <- raw_il8_k[, -1]

# Check that every column name in raw_il8_k (i.e. sample IDs) is found in il8_k_groups$id
if (all(colnames(raw_il8_k) %in% il8_k_groups$id)) {
  cat("All sample names in raw_ are present in _groups$id.\n")
} else {
  missing_cols <- setdiff(colnames(raw_il8_k), il8_k_groups$id)
  cat("The following sample names in raw_ are not found in _groups$id:\n")
  print(missing_cols)
}

# 4. Ensure that the 'group' column is a factor and set the reference level
il8_k_groups$group <- factor(il8_k_groups$group)
il8_k_groups$group <- relevel(il8_k_groups$group, ref = "low")

# 5. Create the DESeqDataSet from the raw counts.
#    Note: We do NOT use tidy=TRUE because our data is in wide format (genes x samples).
dds_il8_k <- DESeqDataSetFromMatrix(countData = raw_il8_k,
                                  colData = il8_k_groups,
                                  design = ~ group)

# Run the DESeq2 analysis (this will estimate size factors and perform normalization internally)
dds_il8_k <- DESeq(dds_il8_k)

# Retrieve the results (differential expression statistics)
res_il8_k <- results(dds_il8_k)
head(as.data.frame(res_il8_k))
summary(res_il8_k)

# Save the results to a CSV file
output_fn <- "C:/Users/brady/OneDrive/Desktop/TCM_RNA/Output/CorrelationTables/il8_k_pairwise.csv"
write.csv(as.data.frame(res), output_fn)

# Optional: Extract normalized counts (if you wish to inspect them later)
norm_counts <- counts(dds, normalized = TRUE)
head(norm_counts)

# Create a volcano plot using adjusted p-values
with(as.data.frame(res_il8_k), plot(log2FoldChange, -log10(padj),
                                  pch = 20, main = "IL-8 (k means)"))

# Highlight significant points consistently using padj
with(subset(as.data.frame(res_il8_k), padj < 0.01 & log2FoldChange >= 2), 
     points(log2FoldChange, -log10(padj), pch = 20, col = "green"))
with(subset(as.data.frame(res_il8_k), padj < 0.01 & log2FoldChange <= -2), 
     points(log2FoldChange, -log10(padj), pch = 20, col = "red"))

#PCA plot
#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation

vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="group") #using the DESEQ2 plotPCA fxn we can