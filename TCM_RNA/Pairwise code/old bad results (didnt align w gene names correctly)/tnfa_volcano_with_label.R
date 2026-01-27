# =========================================
# tnfa High vs. Low: Shrunken Volcano + Highlight
# =========================================

# 0) Libraries ----
library(dplyr)
library(DESeq2)
library(ggplot2)
library(tidyr)
library(tibble)
library(ashr)
library(ggrepel)

# 1) Read data ----
raw         <- read.csv("rawcountdata.csv",    stringsAsFactors = FALSE)
tnfa_groups <- read.csv("tnfa_groups.csv",     stringsAsFactors = FALSE)

# 2) Subset your tnfa samples (adjust cols as needed) ----
raw_tnfa <- raw %>%
  dplyr::select(1,3, 5:7, 9:11, 16:18, 22, 23)

# 3) Move gene IDs into rownames & drop that col ----
rownames(raw_tnfa) <- raw_tnfa[[1]]
raw_tnfa <- raw_tnfa[, -1]

# 4) Sanity‐check sample names ----
if (!all(colnames(raw_tnfa) %in% tnfa_groups$id)) {
  stop("Sample names mismatch!")
} else message("✔ Names OK")

# 5) Prepare metadata ----
tnfa_groups$group <- factor(tnfa_groups$group, levels = c("low","high"))
rownames(tnfa_groups) <- tnfa_groups$id

# 6) Prefilter low‐count genes ----
count_threshold <- 3
min_samples     <- 6
keep <- rowSums(raw_tnfa >= count_threshold) >= min_samples
filtered_raw_tnfa <- raw_tnfa[keep, ]

# 7) DESeq2 + ashr shrinkage ----
dds <- DESeqDataSetFromMatrix(
  countData = filtered_raw_tnfa,
  colData   = tnfa_groups[colnames(filtered_raw_tnfa), ],
  design    = ~ group
)

dds <- DESeq(dds)

res_shrunk <- lfcShrink(
  dds,
  coef = "group_high_vs_low",
  type = "ashr"
)

# 8) Export full results (optional) ----
res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene")
write.csv(res_df, "tnfa_pairwise_filtered.csv", row.names = FALSE)

# 9) Build volcano dataframe ----
volcano_df2 <- res_df %>%
  mutate(
    negLog10padj = -log10(padj),
    direction    = case_when(
      padj < 0.05 & log2FoldChange >=  1 ~ "Up",
      padj < 0.05 & log2FoldChange <= -1 ~ "Down",
      TRUE                               ~ "NS"
    )
  ) %>%
  filter(!is.na(padj), !is.na(log2FoldChange))

# 10) Pick your gene to highlight ----
ens.id   <- "ENSCAFG00000024781"  # ← fill in your ENSCAF ID
symbol   <- "VSIG4"               # ← fill in the label you want
highlight.df <- volcano_df2 %>%
  filter(gene == ens.id) %>%
  mutate(symbol = symbol)

# 11) Draw & save the volcano ----
p_volcano <- ggplot(volcano_df2,
                    aes(x = log2FoldChange, y = negLog10padj)) +
  # main cloud
  geom_point(aes(color = direction,
                 alpha = abs(log2FoldChange)),
             size = 2) +
  scale_color_manual(values = c(NS = "grey70",
                                Down = "red",
                                Up = "green")) +
  scale_alpha(range = c(0.3, 1), guide = "none") +
  # threshold lines
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed",
             color = "darkgrey") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             color = "darkgrey") +
  # circle around your gene
  geom_point(
    data        = highlight.df,
    mapping     = aes(x = log2FoldChange, y = negLog10padj),
    inherit.aes = FALSE,
    shape       = 21,
    fill        = "forestgreen",
    color       = "black",
    size        = 4,
    stroke      = 1
  ) +
  # label it with the human‐readable symbol
  geom_text_repel(
    data        = highlight.df,
    mapping     = aes(x = log2FoldChange,
                      y = negLog10padj,
                      label = symbol),
    inherit.aes = FALSE,
    nudge_x     = 0.2,
    nudge_y     = 0.6,
    size        = 4,
    box.padding = 0.3
  ) +
  # clean theme
  theme_minimal(base_size = 14) +
  labs(
    x = "Shrunken log2 fold change",
    y = "-log10 adjusted p-value"
  ) +
  guides(color = "none")

# plot + save
print(p_volcano)
ggsave("volcano_tnfa_highlight.png",
       p_volcano,
       width  = 7,
       height = 5,
       dpi    = 600,
       bg     = "transparent")
