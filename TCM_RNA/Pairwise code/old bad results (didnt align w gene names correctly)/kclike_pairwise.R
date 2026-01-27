
# =========================================
# kclike High vs. Low Pairwise DEG Analysis
# =========================================

# 0) Load libraries ----
library(dplyr)
library(DESeq2)
library(ggplot2)
library(tidyr)
library(tibble)
library(ashr)

# 1) Read data ----
raw        <- read.csv("rawcountdata.csv", 
                       stringsAsFactors = FALSE)
kclike_groups <- read.csv("kclike_groups.csv", header = TRUE, 
                       stringsAsFactors = FALSE)

# 2) Subset your samples manually (as before) ----
raw_kclike <- raw %>%
  select(1,3:5, 7, 10, 12, 15, 17, 19, 20, 21, 22)

# 3) Move gene IDs into rownames & drop that column ----
rownames(raw_kclike) <- raw_kclike[[1]]
raw_kclike <- raw_kclike[, -1]

# 4) Check sample names against metadata ----
if (!all(colnames(raw_kclike) %in% kclike_groups$id)) {
  stop("Sample names in raw_kclike do not match kclike_groups$id!")
} else {
  message("All sample names match.")
}

# 5) Prepare metadata ----
kclike_groups$group <- factor(kclike_groups$group, levels = c("low","high"))
rownames(kclike_groups) <- kclike_groups$id

# 6) Prefilter low‐count genes ----
count_threshold <- 3
min_samples <- 6  # e.g. half
keep <- rowSums(raw_kclike >= count_threshold) >= min_samples
filtered_raw_kclike <- raw_kclike[keep, ]

# 7) Build & run DESeq2 ----
dds <- DESeqDataSetFromMatrix(
  countData = filtered_raw_kclike,
  colData   = kclike_groups[colnames(filtered_raw_kclike), ],
  design    = ~ group
)

dds <- DESeq(dds)

res_shrunk <- lfcShrink(
  dds,
  coef = "group_high_vs_low",
  type = "ashr"      # <-- swap in ashr
)

# 9) Export results ----
res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene")
write.csv(res_df, "kclike_pairwise_filtered.csv", row.names = FALSE)


# 10) Volcano plot with custom colors & no legend ----
volcano_df <- res_df %>%
  mutate(
    negLog10padj = -log10(padj),
    direction = case_when(
      padj < 0.05 & log2FoldChange >=  1 ~ "Up",
      padj < 0.05 & log2FoldChange <= -1 ~ "Down",
      TRUE                               ~ "NS"
    )
  )

volcano_df2 <- volcano_df %>%
  filter(!is.na(padj), !is.na(log2FoldChange))

p <- ggplot(volcano_df2, aes(x = log2FoldChange, y = negLog10padj)) +
  geom_point(aes(color = direction, alpha = abs(log2FoldChange)), size = 2) +
  scale_color_manual(
    values = c(
      NS   = "grey70",
      Down = "red",
      Up   = "green"
    )
  ) +
  scale_alpha(range = c(0.3, 1), guide = "none") +
  theme_minimal(base_size = 14) +
  labs(
    x     = "Shrunken log2 fold change",
    y     = "-log10 adjusted p-value"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
  guides(color = "none")    # remove the legend

print(p)

ggsave(
  filename = "volcano_kclike.png",
  plot     = p,
  width    = 7,           # in inches
  height   = 5,           # adjust as you like
  dpi      = 600,
  bg       = "transparent"
)
