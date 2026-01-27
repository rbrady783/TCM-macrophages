## =========================================
# Over‐Representation Analysis (ORA) of Top‐5% Genes
# =========================================

# 1) Load libraries ----
library(dplyr)            # data manipulation
library(clusterProfiler)  # ORA functions
library(org.Hs.eg.db)     # human gene annotation
library(ReactomePA)       # Reactome pathway ORA (optional)
library(ggplot2)          # plotting

# 2) Read your differential expression results ----
#    Assumes columns: Hugo_Symbol, log2FoldChange, pvalue, padj, etc.
res <- read.csv("vegf_pairwise_filtered.csv", stringsAsFactors = FALSE)

# 3) Define universe & hit list by top 5% abs(log₂FC) ----
universe_genes <- unique(res$Hugo_Symbol)

n_hits <- ceiling(0.05 * nrow(res))  

sig_genes <- res %>%
  filter(!is.na(log2FoldChange)) %>%        # drop any missing
  arrange(desc(abs(log2FoldChange))) %>%    # rank by magnitude
  slice_head(n = n_hits) %>%                # take top 5%
  pull(Hugo_Symbol)

# 4) Map SYMBOL → ENTREZID (one‐to‐one) ----
univ_map <- bitr(
  universe_genes,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
) %>% distinct(SYMBOL, .keep_all = TRUE)

sig_map <- bitr(
  sig_genes,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
) %>% distinct(SYMBOL, .keep_all = TRUE)

# 5) GO‐Biological Process ORA ----
ego_bp <- enrichGO(
  gene          = sig_map$ENTREZID,
  universe      = univ_map$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.25,
  qvalueCutoff  = 0.25,
  readable      = TRUE
)

# 6) KEGG Pathway ORA ----
ekegg <- enrichKEGG(
  gene          = sig_map$ENTREZID,
  universe      = univ_map$ENTREZID,
  organism      = "hsa",
  keyType       = "kegg",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.25,
  qvalueCutoff  = 0.25
) %>%
  setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# 7) Reactome Pathway ORA (optional) ----
er <- enrichPathway(
  gene          = sig_map$ENTREZID,
  universe      = univ_map$ENTREZID,
  organism      = "human",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.25,
  qvalueCutoff  = 0.25,
  readable      = TRUE
)

# 8) Plot results with checks ----
safe_dotplot <- function(enrich_obj, title) {
  df <- as.data.frame(enrich_obj)
  if (nrow(df) > 0) {
    dotplot(enrich_obj, showCategory = min(10, nrow(df))) +
      ggtitle(title)
  } else {
    message("No enriched terms for: ", title)
  }
}

safe_dotplot(ego_bp,  "GO‐BP Enrichment (top 5 % |log₂FC|)")
safe_dotplot(ekegg,   "KEGG Pathway Enrichment (top 5 % |log₂FC|)")
safe_dotplot(er,      "Reactome Pathway Enrichment (top 5 % |log₂FC|)")

# 9) Save all results to CSV ----
write.csv(as.data.frame(ego_bp),
          "ORA_GO_BP_top5pct.csv",
          row.names = FALSE, quote = FALSE)

write.csv(as.data.frame(ekegg),
          "ORA_KEGG_top5pct.csv",
          row.names = FALSE, quote = FALSE)

write.csv(as.data.frame(er),
          "ORA_Reactome_top5pct.csv",
          row.names = FALSE, quote = FALSE)
