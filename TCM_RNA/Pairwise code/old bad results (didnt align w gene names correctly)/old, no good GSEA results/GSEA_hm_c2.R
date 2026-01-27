library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ggplot2)
library(dplyr)

# 1) Read your preranked DE results ----
pr <- read.csv("xxx_pairwise_filtered.csv", stringsAsFactors = FALSE)

# 2) Clean & build the ranking vector ----
pr2 <- pr %>%
  # drop any rows where the gene or LFC is missing
  filter(!is.na(Hugo_Symbol), !is.na(log2FoldChange)) %>%
  # if you have duplicate Hugo_Symbols, you can collapse them,
  # here we just keep the first occurrence
  distinct(Hugo_Symbol, .keep_all = TRUE)

ranks <- setNames(pr2$log2FoldChange, pr2$Hugo_Symbol) %>%
  sort(decreasing = TRUE)   # now it's guaranteed decreasing

# 3) Grab human Hallmark sets (symbols must match your Hugo_Symbol column) ----
m_df     <- msigdbr(
  species    = "Homo sapiens",
  collection = "H"             # Hallmark
)

term2gene <- m_df %>%
  select(gs_name, gene_symbol)

# 4) Run preranked GSEA ----
gseaRes <- GSEA(
  geneList      = ranks,
  TERM2GENE     = term2gene,
  pvalueCutoff  = 1,       # keep anything with nominal p < 0.25
  pAdjustMethod = "BH",
  verbose       = FALSE
)


# 5) View & plot results ----
head(as.data.frame(gseaRes), 10)

dotplot(gseaRes, showCategory=10) +
  ggtitle("GSEA Hallmark (pre-ranked log₂FC)")

ridgeplot(gseaRes, showCategory=6) +
  ggtitle("Top 6 Hallmarks")

# 6) Save full table ----
write.csv(
  as.data.frame(gseaRes),
  file      = "xxx_GSEA_Hallmark_preRanked.csv",
  row.names = FALSE
)

#
# 4.5) Grab human Canonical (C2) sets ----

m_df_c2     <- msigdbr(
  species    = "Homo sapiens",
  collection = "C2"             # Canonical pathways
)

term2gene_c2 <- m_df_c2 %>%
  select(gs_name, gene_symbol)

# 4.6) Run preranked GSEA on C2 ----
gseaRes_c2 <- GSEA(
  geneList      = ranks,
  TERM2GENE     = term2gene_c2,
  pvalueCutoff  = 1,       # same settings as Hallmark
  pAdjustMethod = "BH",
  verbose       = FALSE
)

# 5.5) View & plot C2 results ----
head(as.data.frame(gseaRes_c2), 10)

dotplot(gseaRes_c2, showCategory=10) +
  ggtitle("GSEA C2 Canonical (pre-ranked log₂FC)")

ridgeplot(gseaRes_c2, showCategory=6) +
  ggtitle("Top 6 C2 Canonical Pathways")

# 6.5) Save full C2 table ----
write.csv(
  as.data.frame(gseaRes_c2),
  file      = "xxx_GSEA_C2_preRanked.csv",
  row.names = FALSE
)
