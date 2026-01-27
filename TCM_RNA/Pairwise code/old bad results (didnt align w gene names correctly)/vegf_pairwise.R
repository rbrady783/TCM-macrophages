
# =========================================
# VEGF High vs. Low Pairwise DEG Analysis
# =========================================

# 0) Load libraries ----
library(dplyr)
library(DESeq2)
library(ggplot2)
library(tidyr)
library(tibble)
library(ashr)

# 1) Read data ----
raw        <- read.csv("rawcountdata.csv", stringsAsFactors = FALSE)
vegf_groups <- read.csv("vegf_groups.csv", header = TRUE, stringsAsFactors = FALSE)

# 2) Subset your samples manually (as before) ----
raw_vegf <- raw %>%
  select(1, 3, 4, 6:8, 14:18, 22, 25)

# 3) Move gene IDs into rownames & drop that column ----
rownames(raw_vegf) <- raw_vegf[[1]]
raw_vegf <- raw_vegf[, -1]

# 4) Check sample names against metadata ----
if (!all(colnames(raw_vegf) %in% vegf_groups$id)) {
  stop("Sample names in raw_vegf do not match vegf_groups$id!")
} else {
  message("All sample names match.")
}

# 5) Prepare metadata ----
vegf_groups$group <- factor(vegf_groups$group, levels = c("low","high"))
rownames(vegf_groups) <- vegf_groups$id

# 6) Prefilter low‐count genes ----
count_threshold <- 3
min_samples <- 6  # e.g. half
keep <- rowSums(raw_vegf >= count_threshold) >= min_samples
filtered_raw_vegf <- raw_vegf[keep, ]

# 7) Build & run DESeq2 ----
dds <- DESeqDataSetFromMatrix(
  countData = filtered_raw_vegf,
  colData   = vegf_groups[colnames(filtered_raw_vegf), ],
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
write.csv(res_df, "vegf_pairwise_filtered.csv", row.names = FALSE)


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
  filename = "volcano_VEGF.png",
  plot     = p,
  width    = 7,           # in inches
  height   = 5,           # adjust as you like
  dpi      = 600,
  bg       = "transparent"
)
