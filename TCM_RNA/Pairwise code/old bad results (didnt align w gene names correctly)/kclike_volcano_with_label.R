# =========================================
# kclike High vs. Low: Shrunken Volcano + Highlight
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
kclike_groups <- read.csv("kclike_groups.csv",     stringsAsFactors = FALSE)

# 2) Subset your kclike samples (adjust cols as needed) ----
raw_kclike <- raw %>%
  dplyr::select(1,3:5, 7, 10, 12, 15, 17, 19, 20, 21, 22)

# 3) Move gene IDs into rownames & drop that col ----
rownames(raw_kclike) <- raw_kclike[[1]]
raw_kclike <- raw_kclike[, -1]

# 4) Sanity‐check sample names ----
if (!all(colnames(raw_kclike) %in% kclike_groups$id)) {
  stop("Sample names mismatch!")
} else message("✔ Names OK")

# 5) Prepare metadata ----
kclike_groups$group <- factor(kclike_groups$group, levels = c("low","high"))
rownames(kclike_groups) <- kclike_groups$id

# 6) Prefilter low‐count genes ----
count_threshold <- 3
min_samples     <- 6
keep <- rowSums(raw_kclike >= count_threshold) >= min_samples
filtered_raw_kclike <- raw_kclike[keep, ]

# 7) DESeq2 + ashr shrinkage ----
dds <- DESeqDataSetFromMatrix(
  countData = filtered_raw_kclike,
  colData   = kclike_groups[colnames(filtered_raw_kclike), ],
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
write.csv(res_df, "kclike_pairwise_filtered.csv", row.names = FALSE)

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
ens.id   <- "ENSCAFG00000016041"  # ← fill in your ENSCAF ID
symbol   <- "CRYBG2"               # ← fill in the label you want
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
    nudge_x     = -0.3,
    nudge_y     = 0.3,
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
ggsave("volcano_kclike_highlight.png",
       p_volcano,
       width  = 7,
       height = 5,
       dpi    = 600,
       bg     = "transparent")
