# 1) Read your preranked DE results ----
pr <- read.csv("ccl2_pairwise_filtered.csv", stringsAsFactors = FALSE)

# 2) Clean & build the ranking vector ----
#library(dplyr)
pr2 <- pr %>%
  filter(!is.na(Hugo_Symbol), !is.na(log2FoldChange)) %>%
  distinct(Hugo_Symbol, .keep_all = TRUE)

ranks <- setNames(pr2$log2FoldChange, pr2$Hugo_Symbol) %>%
  sort(decreasing = TRUE)

# 3) Pull the C5 (GO) gene‐sets ----
#library(msigdbr)
#m_df_c5 <- msigdbr(
 # species    = "Homo sapiens",
  #collection = "C5"
#)

#term2gene_c5 <- m_df_c5 %>%
 # dplyr::select(gs_name, gene_symbol)

# 4) Run preranked GSEA on C5 ----
#library(clusterProfiler)

gsea_c5 <- GSEA(
  geneList      = ranks,
  TERM2GENE     = term2gene_c5,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  verbose       = FALSE
)

res_c5_sig <- as.data.frame(gsea_c5)

# 5) Plot top GO hits ----
dotplot(gsea_c5, showCategory = head(res_c5_sig$ID, 8)) +
  ggtitle("GSEA: Top C5 (GO) Terms")

# 6) Save results ----
#out_dir  <- "C:/Users/brady/OneDrive/Desktop/TCM_RNA/Pairwise code/GSEA results"
out_file <- file.path(out_dir, "ccl2_GSEA_C5_preRanked.csv")

write.csv(
  res_c5_sig,
  file      = out_file,
  row.names = FALSE
)
