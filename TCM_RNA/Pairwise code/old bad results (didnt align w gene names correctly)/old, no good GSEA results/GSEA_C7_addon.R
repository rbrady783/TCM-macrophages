

# 1) Read your preranked DE results ----
pr <- read.csv("tgfb_pairwise_filtered.csv", stringsAsFactors = FALSE)

# 2) Clean & build the ranking vector ----
pr2 <- pr %>%
  # drop any rows where the gene or LFC is missing
  filter(!is.na(Hugo_Symbol), !is.na(log2FoldChange)) %>%
  # if you have duplicate Hugo_Symbols, you can collapse them,
  # here we just keep the first occurrence
  distinct(Hugo_Symbol, .keep_all = TRUE)

ranks <- setNames(pr2$log2FoldChange, pr2$Hugo_Symbol) %>%
  sort(decreasing = TRUE)   # now it's guaranteed decreasing


# 1) Pull the C7 immunologic signatures
#m_df_c7 <- msigdbr(
  #species    = "Homo sapiens",
  #collection = "C7"
#) 

#term2gene_c7 <- m_df_c7 %>% 
  #select(gs_name, gene_symbol)

# 2) Run preranked GSEA on C7
gsea_c7 <- GSEA(
  geneList      = ranks,
  TERM2GENE     = term2gene_c7,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  verbose       = FALSE
)

res_c7_sig <- as.data.frame(gsea_c7) 

# 4) Plot top immune‐related hits
dotplot(gsea_c7, showCategory = head(res_c7_sig$ID, 8))

write.csv(gsea_c7,
  file      = "tgfb_GSEA_C7_preRanked.csv",
  row.names = FALSE
)

