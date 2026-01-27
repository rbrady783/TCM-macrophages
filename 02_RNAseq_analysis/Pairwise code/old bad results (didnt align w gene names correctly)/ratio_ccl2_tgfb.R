# =========================================
# ratio High vs. Low Pairwise DEG Analysis
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
ratio_groups <- read.csv("ratio_groups.csv", header = TRUE, 
                        stringsAsFactors = FALSE)

# 2) Subset your samples manually (as before) ----
raw_ratio <- raw %>%
  select(1,3,5:8,10,13,14,16,19,20,25)

# 3) Move gene IDs into rownames & drop that column ----
rownames(raw_ratio) <- raw_ratio[[1]]
raw_ratio <- raw_ratio[, -1]

# 4) Check sample names against metadata ----
if (!all(colnames(raw_ratio) %in% ratio_groups$id)) {
  stop("Sample names in raw_ratio do not match ratio_groups$id!")
} else {
  message("All sample names match.")
}

# 5) Prepare metadata ----
ratio_groups$group <- factor(ratio_groups$group, levels = c("low","high"))
rownames(ratio_groups) <- ratio_groups$id

# 6) Prefilter low‐count genes ----
count_threshold <- 3
min_samples <- 6  # e.g. half
keep <- rowSums(raw_ratio >= count_threshold) >= min_samples
filtered_raw_ratio <- raw_ratio[keep, ]

# 7) Build & run DESeq2 ----
dds <- DESeqDataSetFromMatrix(
  countData = filtered_raw_ratio,
  colData   = ratio_groups[colnames(filtered_raw_ratio), ],
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
write.csv(res_df, "ratio_pairwise_filtered.csv", row.names = FALSE)


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

library(ggplot2)
library(ggrepel)    # for nice text labels

# 1) Decide your cap
ycap <- 6

# 2) Prepare your plotting frame
volcano_df2 <- volcano_df2 %>%
  mutate(
    above_cap = negLog10padj > ycap,
    y_plot    = ifelse(above_cap, ycap, negLog10padj)
  )

# 3) Find the single outlier row
outlier <- volcano_df2 %>% filter(above_cap) %>% slice_max(negLog10padj, n = 1)

# 4) Plot without arrow, with asterisk annotation
p_volcano <- ggplot(volcano_df2, aes(x = log2FoldChange, y = y_plot)) +
  geom_point(aes(color = direction), size = 2, alpha = 0.8) +
  # dashed thresholds
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
  # asterisk for the single capped outlier
  geom_text(
    data    = outlier,
    aes(label = "*"),
    vjust   = -0.1,
    size    = 6,
    color   = "black"
  ) +
  scale_color_manual(values = c(NS = "grey70", Down = "red", Up = "green")) +
  coord_cartesian(ylim = c(0, ycap + 1)) +
  theme_minimal(base_size = 14) +
  labs(
    x     = "Shrunken log2 fold change",
    y     = "-log10 adjusted p-value"
  ) +
  guides(color = "none")

# 5) Display
print(p_volcano)

ggsave(
  filename = "volcano_ratio.png",
  plot     = p_volcano,
  width    = 7,           # in inches
  height   = 5,           # adjust as you like
  dpi      = 600,
  bg       = "transparent"
)