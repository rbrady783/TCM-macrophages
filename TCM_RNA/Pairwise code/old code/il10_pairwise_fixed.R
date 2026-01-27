# =========================================
# il10 High vs. Low Pairwise DEG Analysis
# =========================================

# 0) Load libraries ----
library(dplyr)
library(DESeq2)
library(ggplot2)
library(tibble)
library(ashr)

# 1) Read data ----
raw        <- read.csv("rawcountdata.csv", 
                       stringsAsFactors = FALSE)
il10_groups <- read.csv("il10_groups.csv", header = TRUE, 
                       stringsAsFactors = FALSE)

# 2) Subset your samples manually (as before) ----
raw_il10 <- raw %>%
  select(1,7:9,11, 14:18, 22, 23, 25)

# 3) Move gene IDs into rownames & drop that column ----
rownames(raw_il10) <- raw_il10[[1]]
raw_il10 <- raw_il10[, -1]

# 4) Check sample names against metadata ----
if (!all(colnames(raw_il10) %in% il10_groups$id)) {
  stop("Sample names in raw_il10 do not match il10_groups$id!")
} else {
  message("All sample names match.")
}

# 5) Prepare metadata ----
il10_groups$group <- factor(il10_groups$group, levels = c("low","high"))
rownames(il10_groups) <- il10_groups$id

# 6) Prefilter low‐count genes ----
count_threshold <- 3
min_samples <- 6  # e.g. half
keep <- rowSums(raw_il10 >= count_threshold) >= min_samples
filtered_raw_il10 <- raw_il10[keep, ]

# 7) Build & run DESeq2 ----
dds_il10 <- DESeqDataSetFromMatrix(
  countData = filtered_raw_il10,
  colData   = il10_groups[colnames(filtered_raw_il10), ],
  design    = ~ group
)

dds_il10 <- DESeq(dds_il10)

res_shrunk_il10 <- lfcShrink(
  dds_il10,
  coef = "group_high_vs_low",
  type = "ashr"      # <-- swap in ashr
)

# 8) Export results ----
res_df_il10 <- as.data.frame(res_shrunk_il10) %>%
  rownames_to_column("gene")
write.csv(res_df_il10, "il10_pairwise_fixed.csv", row.names = FALSE)

# 9) Build volcano dataframe ----
volcano_df2 <- res_df_il10 %>%
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
ens.id   <- "ENSCAFG00000018219"  # ← fill in your ENSCAF ID
symbol   <- "CSF1R"               # ← fill in the label you want
highlight.df <- volcano_df2 %>%
  filter(gene == ens.id) %>%
  mutate(symbol = symbol)

# 11) Draw & save the volcano (bolder, no axis labels, corner tag) ----

panel_tag <- "IL-10"  # <- change to "CCL2", "IL-8", etc. for other panels

p_volcano <- ggplot(
  volcano_df2,
  aes(x = log2FoldChange, y = negLog10padj)
) +
  # main cloud: bigger points, high alpha (bolder)
  geom_point(aes(color = direction),
             size = 2.6, alpha = 0.9) +
  
  # high-contrast, colorblind-friendly palette
  scale_color_manual(values = c(NS = "#9E9E9E",
                                Down = "red",
                                Up   = "green")) +
  
  # threshold lines: darker + slightly thicker
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed",
             color = "grey30",
             linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             color = "grey30",
             linewidth = 0.5) +
  
  # panel tag in upper-left (inside plot area)
  annotate("text", x = -Inf, y = Inf, label = panel_tag,
           hjust = -0.1, vjust = 1.2, size = 5, fontface = "bold") +
  
  # highlight your gene of interest (same as before)
  geom_point(
    data        = highlight.df,
    mapping     = aes(x = log2FoldChange, y = negLog10padj),
    inherit.aes = FALSE,
    shape       = 21,
    fill        = "gold",
    color       = "black",
    size        = 4,
    stroke      = 1
  ) +
  ggrepel::geom_text_repel(
    data        = highlight.df,
    mapping     = aes(x = log2FoldChange, y = negLog10padj, label = symbol),
    inherit.aes = FALSE,
    nudge_x     = 0.2,
    nudge_y     = 0.4,
    size        = 4,
    box.padding = 0.3
  ) +
  
  # crisp figure style: strong axes & border, no axis titles
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.title   = element_blank()
  ) +
  guides(color = "none")

# plot + save
print(p_volcano)

ggsave("volcano_IL10.pdf", p_volcano, width = 7, height = 5, device = cairo_pdf)
ggsave("volcano_IL10.png", p_volcano, width = 7, height = 5, dpi = 600)
