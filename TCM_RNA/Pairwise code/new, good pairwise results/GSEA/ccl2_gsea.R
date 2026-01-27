library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ggplot2)
library(dplyr)

# 1) Read your preranked DE results ----change cytokine as needed----
pr <- read.csv("ccl2_quartile_results.csv", stringsAsFactors = FALSE) 

# 2) Clean & build the ranking vector ----
pr2 <- pr %>%
  # drop any rows where the gene or LFC is missing
  filter(!is.na(Hugo_Symbol), !is.na(log2FoldChange)) %>%
  # if you have duplicate Hugo_Symbols, you can collapse them,
  # here we just keep the first occurrence
  distinct(Hugo_Symbol, .keep_all = TRUE)

ranks <- setNames(pr2$log2FoldChange, pr2$Hugo_Symbol) %>%
  sort(decreasing = TRUE)   # now it's guaranteed decreasing

# Quality control check
cat("Total genes in ranking:", length(ranks), "\n")
hist(ranks, breaks=50, main="Distribution of log2FC rankings")

# 3) Grab human Hallmark sets (H) ----
m_df_h <- msigdbr(
  species    = "Homo sapiens",
  collection = "H"             # Hallmark
)
term2gene_h <- m_df_h %>%
  select(gs_name, gene_symbol)

# Check gene overlap
cat("Genes in Hallmark sets:", length(unique(m_df_h$gene_symbol)), "\n")
cat("Overlap with ranking:", length(intersect(names(ranks), unique(m_df_h$gene_symbol))), "\n")

# Set seed before running GSEA
set.seed(123)  # or any number you prefer

gseaRes_h <- GSEA(
  geneList      = ranks,
  TERM2GENE     = term2gene_h,
  pvalueCutoff  = 0.25,
  pAdjustMethod = "BH",
  minGSSize     = 10,
  maxGSSize     = 1000,
  eps           = 0,
  seed          = TRUE,    # Now this should work
  verbose       = FALSE
)

# 5) View & plot Hallmark results ----
head(as.data.frame(gseaRes_h), 10)
dotplot(gseaRes_h, showCategory=10) +
  ggtitle("GSEA Hallmark (pre-ranked log₂FC)")
ridgeplot(gseaRes_h, showCategory=6) +
  ggtitle("Top 6 Hallmarks")

# 6) Save Hallmark table ----
write.csv(
  as.data.frame(gseaRes_h),
  file      = "ccl2_GSEA_Hallmark_preRanked.csv",
  row.names = FALSE
)

# Look for strong effect sizes
sig_results <- as.data.frame(gseaRes_h) %>%
  filter(p.adjust < 0.25) %>%
  arrange(desc(abs(NES)))

head(sig_results[,c("Description", "NES", "pvalue", "p.adjust")])

# 7) Grab human C2 Canonical sets ----
m_df_c2 <- msigdbr(
  species    = "Homo sapiens",
  collection = "C2"             # Canonical pathways
)
term2gene_c2 <- m_df_c2 %>%
  select(gs_name, gene_symbol)

# 8) Run preranked GSEA - C2 ----
gseaRes_c2 <- GSEA(
  geneList      = ranks,
  TERM2GENE     = term2gene_c2,
  pvalueCutoff  = 0.25,       # same settings as Hallmark
  pAdjustMethod = "BH",
  minGSSize     = 10,
  maxGSSize     = 1000,
  eps           = 0,
  seed          = TRUE,
  verbose       = FALSE
)

# 9) View & plot C2 results ----
head(as.data.frame(gseaRes_c2), 10)
dotplot(gseaRes_c2, showCategory=10) +
  ggtitle("GSEA C2 Canonical (pre-ranked log₂FC)")
ridgeplot(gseaRes_c2, showCategory=6) +
  ggtitle("Top 6 C2 Canonical Pathways")

# 10) Save C2 table ----
write.csv(
  as.data.frame(gseaRes_c2),
  file      = "ccl2_GSEA_C2_preRanked.csv",
  row.names = FALSE
)


# 15) Grab human C6 Oncogenic sets ----
m_df_c6 <- msigdbr(
  species    = "Homo sapiens",
  collection = "C6"             # Oncogenic signatures
)
term2gene_c6 <- m_df_c6 %>%
  select(gs_name, gene_symbol)

# 16) Run preranked GSEA - C6 ----
gseaRes_c6 <- GSEA(
  geneList      = ranks,
  TERM2GENE     = term2gene_c6,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  minGSSize     = 10,
  maxGSSize     = 1000,
  eps           = 0,
  seed          = TRUE,
  verbose       = FALSE
)

# 17) View & plot C6 results ----
head(as.data.frame(gseaRes_c6), 10)
dotplot(gseaRes_c6, showCategory=10) +
  ggtitle("GSEA C6 Oncogenic (pre-ranked log₂FC)")
ridgeplot(gseaRes_c6, showCategory=6) +
  ggtitle("Top 6 Oncogenic Signatures")

# 18) Save C6 table ----
write.csv(
  as.data.frame(gseaRes_c6),
  file      = "ccl2_GSEA_C6_preRanked.csv",
  row.names = FALSE
)

# 19) Grab human C7 Immunologic sets ----
m_df_c7 <- msigdbr(
  species    = "Homo sapiens",
  collection = "C7"             # Immunologic signatures
)
term2gene_c7 <- m_df_c7 %>%
  select(gs_name, gene_symbol)

# 20) Run preranked GSEA - C7 ----
gseaRes_c7 <- GSEA(
  geneList      = ranks,
  TERM2GENE     = term2gene_c7,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  minGSSize     = 10,
  maxGSSize     = 1000,
  eps           = 0,
  seed          = TRUE,
  verbose       = FALSE
)

# 21) View & plot C7 results ----
head(as.data.frame(gseaRes_c7), 10)
dotplot(gseaRes_c7, showCategory=10) +
  ggtitle("GSEA C7 Immunologic (pre-ranked log₂FC)")
ridgeplot(gseaRes_c7, showCategory=6) +
  ggtitle("Top 6 Immunologic Signatures")

# 22) Save C7 table ----
write.csv(
  as.data.frame(gseaRes_c7),
  file      = "ccl2_GSEA_C7_preRanked.csv",
  row.names = FALSE
)

# 23) Combine all results into master table ----
# Extract cytokine name from filename (adjust as needed)
cytokine_name <- "CCL2"  # Change this to your actual cytokine name

collections <- c("Hallmark", "C2", "C6", "C7")
gsea_list <- list(gseaRes_h, gseaRes_c2, gseaRes_c6, gseaRes_c7)
names(gsea_list) <- collections

# Create master results table
master_results <- data.frame()

for(i in 1:length(collections)) {
  # Get results for this collection
  collection_results <- as.data.frame(gsea_list[[i]])
  
  # Add identifying columns
  collection_results$Cytokine <- cytokine_name
  collection_results$Collection <- collections[i]
  
  # Reorder columns to put identifiers first
  collection_results <- collection_results[, c("Cytokine", "Collection", 
                                               names(collection_results)[!names(collection_results) %in% c("Cytokine", "Collection")])]
  
  # Combine with master table
  master_results <- rbind(master_results, collection_results)
}

# Sort by significance
master_results <- master_results[order(master_results$p.adjust), ]

print(paste("Total pathways analyzed:", nrow(master_results)))
print(paste("Significant at FDR < 0.25:", sum(master_results$p.adjust < 0.25, na.rm = TRUE)))

# 24) Save master results table ----
write.csv(
  master_results,
  file      = paste0(cytokine_name, "_GSEA_All_Collections.csv"),
  row.names = FALSE
)
