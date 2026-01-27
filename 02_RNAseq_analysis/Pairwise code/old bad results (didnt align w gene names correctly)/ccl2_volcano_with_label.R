# ==============================================
# CCL2 High vs. Low: DEG → Shrunk Volcano w/ Highlight
# ==============================================

# 0) load all required packages ----
library(dplyr)
library(DESeq2)
library(ggplot2)
library(tidyr)
library(tibble)
library(ashr)
library(ggrepel)

# 1) read in counts & metadata ----
raw        <- read.csv("rawcountdata.csv",   stringsAsFactors = FALSE)
ccl2_groups <- read.csv("ccl2_groups.csv",   stringsAsFactors = FALSE)

# 2) subset to your CCL2 samples ----
raw_ccl2 <- raw %>%
  dplyr::select(1,3,4,6,11,13,15,16,18:20,22,25)

# 3) move gene IDs into rownames & drop first col ----
rownames(raw_ccl2) <- raw_ccl2[[1]]
raw_ccl2 <- raw_ccl2[,-1]

# 4) sanity‐check sample names ----
if (!all(colnames(raw_ccl2) %in% ccl2_groups$id)) {
  stop("Sample names in raw_ccl2 do not match ccl2_groups$id!")
} else {
  message("✔ Sample names match.")
}

# 5) prep the colData for DESeq2 ----
ccl2_groups$group <- factor(ccl2_groups$group, levels = c("low","high"))
rownames(ccl2_groups) <- ccl2_groups$id

# 6) prefilter low‐count genes ----
count_threshold <- 3
min_samples     <- 6
keep <- rowSums(raw_ccl2 >= count_threshold) >= min_samples
filtered_raw_ccl2 <- raw_ccl2[keep, ]

# 7) run DESeq2 + ashr shrinkage ----
dds <- DESeqDataSetFromMatrix(
  countData = filtered_raw_ccl2,
  colData   = ccl2_groups[colnames(filtered_raw_ccl2), ],
  design    = ~ group
)
dds <- DESeq(dds)
res_shrunk <- lfcShrink(
  dds,
  coef = "group_high_vs_low",
  type = "ashr"
)

# 8) export full results if you like ----
res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene")
write.csv(res_df,
          "ccl2_pairwise_filtered.csv",
          row.names = FALSE)

# 9) build volcano data.frame ----
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

# 10) cap any extreme –log10(padj) outlier at ycap, mark with asterisk ----
ycap <- 6
volcano_df2 <- volcano_df2 %>%
  mutate(
    above_cap = negLog10padj > ycap,
    y_plot    = ifelse(above_cap, ycap, negLog10padj)
  )
outlier <- volcano_df2 %>%
  filter(above_cap) %>%
  slice_max(negLog10padj, n = 1)


# 11) pick your ENSCAF ID and give it the human‐readable name
ens.id   <- "ENSCAFG00000015233"
readable <- "EWSR1"
highlight.df <- volcano_df2 %>%
  filter(gene == ens.id) %>%
  mutate(symbol = readable)

highlight.df

# 12) now draw the volcano, layering on the asterisk + circle + label
p_volcano <- ggplot(volcano_df2, aes(x = log2FoldChange, y = y_plot)) +
  # main points
  geom_point(aes(color = direction, alpha = abs(log2FoldChange)), size = 2) +
  
  # dashed thresholds
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
  
  # asterisk for the capped outlier
  geom_text(
    data        = outlier,
    mapping     = aes(x = log2FoldChange, y = y_plot, label = "*"),
    inherit.aes = FALSE,
    vjust       = -0.2,
    size        = 6
  ) +
  
  # yellow circle around the ENSCAF point
  geom_point(
    data        = highlight.df,
    mapping     = aes(x = log2FoldChange, y = y_plot),
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
    mapping     = aes(x = log2FoldChange, y = y_plot, label = symbol),
    inherit.aes = FALSE,
    nudge_x     = 0.5,   # move right by 0.5 on the x‐axis
    nudge_y     = 0.7,   # move up by 1 on the y‐axis
    size        = 4,
    box.padding = 0.3
  ) +
  
  # finishes
  scale_color_manual(values = c(NS = "grey70", Down = "red", Up = "green")) +
  scale_alpha(range = c(0.3,1), guide="none") +
  coord_cartesian(ylim = c(0, ycap + 1)) +
  theme_minimal(base_size = 14) +
  labs(y="-log10 adjusted p-value") +
  guides(color="none")

print(p_volcano)

ggsave("volcano_ccl2_EWSR1.png", p_volcano, width=7, height=5, dpi=600)
