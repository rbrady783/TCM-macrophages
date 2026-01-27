#Pairwise comparison
###########################
#library(htmltools)
#library(ggplot2)
#library(DESeq2)
#library(dplyr)

# Download data
raw <- read.csv('RawCountData_CellLines_STAR_GeneCounts_RNAseq_For_Rachel.csv')
head(raw)

# Select relevant columns
raw_il8_low <- select(raw, 1,3,5,9,16,18:23,25)

# Groups
il8_groups_low <- read.csv('il8_groups_low.cutoff.csv')

# Construct DESEQDataSet Object
dds <- DESeqDataSetFromMatrix(countData=raw_il8_low, 
                              colData=il8_groups_low, 
                              design=~group, tidy = TRUE)

# run DEQSEQ
dds <- DESeq(dds)

# look at results
res <- results(dds)
head(results(dds, tidy=TRUE))

summary(res)

# save results
fn <- paste0("C:/Users/brady/OneDrive/Desktop/TCM_RNA/Output/CorrelationTables/il8_pairwise_low.cutoff.csv")
write.csv(res, fn)

# plots
plotCounts(dds, gene="ENSCAFG00000014882", intgroup="group")

#volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="IL-8 Lower cutoff", ylim=c(0,10), xlim=c(-10,10)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<=.01 & log2FoldChange >=2), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<=.01 & log2FoldChange <=-2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#PCA plot
#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation

vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="group") #using the DESEQ2 plotPCA fxn we can
