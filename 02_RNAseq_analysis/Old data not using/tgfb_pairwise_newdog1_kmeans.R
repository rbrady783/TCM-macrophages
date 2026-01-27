# Select relevant columns
raw_tgfb_k <- select(raw, 1,3,4,6,7,13,14,19,21:23)

# 3. Read the sample metadata from tgfb_k_groups.csv
tgfb_k_groups <- read.csv('tgfb_k_groups.csv', header = TRUE, stringsAsFactors = FALSE)

# Set the row names to the gene IDs (first column), then remove that column.
rownames(raw_tgfb_k) <- raw_tgfb_k[[1]]
raw_tgfb_k <- raw_tgfb_k[, -1]

# Check that every column name in raw_tgfb_k (i.e. sample IDs) is found in tgfb_k_groups$id
if (all(colnames(raw_tgfb_k) %in% tgfb_k_groups$id)) {
  cat("All sample names in raw_ are present in _groups$id.\n")
} else {
  missing_cols <- setdiff(colnames(raw_tgfb_k), tgfb_k_groups$id)
  cat("The following sample names in raw_ are not found in _groups$id:\n")
  print(missing_cols)
}

# 4. Ensure that the 'group' column is a factor and set the reference level
tgfb_k_groups$group <- factor(tgfb_k_groups$group)
tgfb_k_groups$group <- relevel(tgfb_k_groups$group, ref = "low")

# 5. Create the DESeqDataSet from the raw counts.
#    Note: We do NOT use tidy=TRUE because our data is in wide format (genes x samples).
dds_tgfb_k <- DESeqDataSetFromMatrix(countData = raw_tgfb_k,
                                   colData = tgfb_k_groups,
                                   design = ~ group)

# Run the DESeq2 analysis (this will estimate size factors and perform normalization internally)
dds_tgfb_k <- DESeq(dds_tgfb_k)

# Retrieve the results (differential expression statistics)
res_tgfb_k <- results(dds_tgfb_k)
head(as.data.frame(res_tgfb_k))
summary(res_tgfb_k)

# Save the results to a CSV file
output_fn <- "C:/Users/brady/OneDrive/Desktop/TCM_RNA/Output/CorrelationTables/tgfb_k_pairwise.csv"
write.csv(as.data.frame(res), output_fn)

# Optional: Extract normalized counts (if you wish to inspect them later)
norm_counts <- counts(dds, normalized = TRUE)
head(norm_counts)

# Create a volcano plot using adjusted p-values
with(as.data.frame(res_tgfb_k), plot(log2FoldChange, -log10(padj),
                                   pch = 20, main = "TGFb (k means)"))

# Highlight significant points consistently using padj
with(subset(as.data.frame(res_tgfb_k), padj < 0.01 & log2FoldChange >= 2), 
     points(log2FoldChange, -log10(padj), pch = 20, col = "green"))
with(subset(as.data.frame(res_tgfb_k), padj < 0.01 & log2FoldChange <= -2), 
     points(log2FoldChange, -log10(padj), pch = 20, col = "red"))

#PCA plot
#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation

vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="group") #using the DESEQ2 plotPCA fxn we can