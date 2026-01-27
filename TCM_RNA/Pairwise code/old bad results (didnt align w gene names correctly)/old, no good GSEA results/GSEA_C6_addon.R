# ─────────────────────────────────────────────────────────────────────────────
# Preranked GSEA against MSigDB C6 (Oncogenic Signatures)
# ─────────────────────────────────────────────────────────────────────────────

# (Make sure you’ve already done:)
# library(dplyr)
# library(clusterProfiler)
# library(msigdbr)
# library(ggplot2)

# 1) Read your preranked DE results ----
pr <- read.csv("ccl2_pairwise_filtered.csv", stringsAsFactors = FALSE)

# 2) Clean & build the ranking vector ----
pr2 <- pr %>%
  # drop any rows where the gene or LFC is missing
  filter(!is.na(Hugo_Symbol), !is.na(log2FoldChange)) %>%
  # if you have duplicate Hugo_Symbols, you can collapse them,
  # here we just keep the first occurrence
  distinct(Hugo_Symbol, .keep_all = TRUE)

ranks <- setNames(pr2$log2FoldChange, pr2$Hugo_Symbol) %>%
  sort(decreasing = TRUE)   # now it's guaranteed decreasing

# 1) pull down the C6 gene‐sets
#m_df_c6 <- msigdbr(
 # species    = "Homo sapiens",
  #collection = "C6"
#)

#term2gene_c6 <- m_df_c6 %>%
 #dplyr::select(
  #gs_name,      # the pathway ID / name
  #gene_symbol  # the HGNC symbol
  #)

# 2) run preranked GSEA
gsea_c6 <- GSEA(
  geneList      = ranks,          # from your preranked log₂FC vector
  TERM2GENE     = term2gene_c6,
  pvalueCutoff  = 1,              # keep everything, filter later
  pAdjustMethod = "BH",
  verbose       = FALSE
)

# 3) pull out the results
res_c6_sig <- as.data.frame(gsea_c6)

# 6) Save results ----
out_dir  <- "C:/Users/brady/OneDrive/Desktop/TCM_RNA/Pairwise code/GSEA results"
out_file <- file.path(out_dir, "ccl2_GSEA_C6_preRanked.csv")

write.csv(
  res_c6_sig,
  file      = out_file,
  row.names = FALSE
)

# 5) quick dot‐plot of the top hits
#    (you can tweak showCategory = N to show more/fewer)
dotplot(
  gsea_c6,
  showCategory = head(res_c6_sig$ID, 8)
) 
