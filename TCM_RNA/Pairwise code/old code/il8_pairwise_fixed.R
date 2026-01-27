# =========================================
# il8 High vs. Low Pairwise DEG Analysis
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
il8_groups <- read.csv("il8_groups.csv", header = TRUE, 
                        stringsAsFactors = FALSE)

# 2) Subset your samples manually (as before) ----
raw_il8 <- raw %>%
  select(1,3,5,9,16:23,25)

# 3) Move gene IDs into rownames & drop that column ----
rownames(raw_il8) <- raw_il8[[1]]
raw_il8 <- raw_il8[, -1]

# 4) Check sample names against metadata ----
if (!all(colnames(raw_il8) %in% il8_groups$id)) {
  stop("Sample names in raw_il8 do not match il8_groups$id!")
} else {
  message("All sample names match.")
}

# 5) Prepare metadata ----
il8_groups$group <- factor(il8_groups$group, levels = c("low","high"))
rownames(il8_groups) <- il8_groups$id

# 6) Prefilter low‐count genes ----
count_threshold <- 3
min_samples <- 6  # e.g. half
keep <- rowSums(raw_il8 >= count_threshold) >= min_samples
filtered_raw_il8 <- raw_il8[keep, ]

# 7) Build & run DESeq2 ----
dds_il8 <- DESeqDataSetFromMatrix(
  countData = filtered_raw_il8,
  colData   = il8_groups[colnames(filtered_raw_il8), ],
  design    = ~ group
)

dds_il8 <- DESeq(dds_il8)

res_shrunk_il8 <- lfcShrink(
  dds_il8,
  coef = "group_high_vs_low",
  type = "ashr"      # <-- swap in ashr
)

# 8) Export results ----
res_df_il8 <- as.data.frame(res_shrunk_il8) %>%
  rownames_to_column("gene")
write.csv(res_df_il8, "il8_pairwise_fixed.csv", row.names = FALSE)

# 9) Build a volcano data.frame (no capping/outlier handling) ------------------
# - Compute -log10(padj)
# - Classify genes as Up / Down / NS at padj<0.05 and |log2FC|≥1
# - Only NS points are semi-transparent to reduce visual fog
volcano_df2 <- res_df %>%
  dplyr::mutate(
    negLog10padj = -log10(padj),
    direction    = dplyr::case_when(
      padj < 0.05 & log2FoldChange >=  1 ~ "Up",
      padj < 0.05 & log2FoldChange <= -1 ~ "Down",
      TRUE                               ~ "NS"
    ),
    alpha_pt     = ifelse(direction == "NS", 0.25, 1)  # NS points lighter
  ) %>%
  dplyr::filter(!is.na(padj), !is.na(log2FoldChange))

# Option: a little headroom for the y-axis
y_max <- max(volcano_df2$negLog10padj, na.rm = TRUE)

# 10) (Optional) highlight a gene of interest ----------------------------------
# Replace 'ens.id' and 'readable' as needed; leave as-is to skip highlighting
ens.id   <- "ENSCAFG00000018348"  # example: SNAI1
readable <- "CCL7"
highlight.df <- volcano_df2 %>%
  dplyr::filter(gene == ens.id) %>%
  dplyr::mutate(symbol = readable)

# 11) Draw the volcano (no outlier marker, no axis titles) ---------------------
panel_tag <- "IL-8"   # change to "CCL2" etc. if desired

p_volcano <- ggplot(volcano_df2, aes(x = log2FoldChange, y = negLog10padj)) +
  # main scatter; only NS is semi-transparent
  geom_point(aes(color = direction, alpha = alpha_pt),
             size = 2.2, show.legend = FALSE) +
  
  # dashed thresholds (|log2FC|=1 and padj=0.05)
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed", color = "grey30", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "grey30", linewidth = 0.4) +
  
  # panel tag in the upper-left corner
  annotate("text", x = -Inf, y = Inf, label = panel_tag,
           hjust = -0.1, vjust = 1.2, size = 5, fontface = "bold") +
  
  # optional highlight point + label (draws nothing if highlight.df is empty)
  geom_point(
    data  = highlight.df,
    shape = 21, fill = "gold", color = "black",
    size  = 3.8, stroke = 0.6
  ) +
  ggrepel::geom_text_repel(
    data = highlight.df, aes(label = symbol),
    nudge_x = 0.6, nudge_y = 0.8, size = 4,
    segment.color = "black", segment.size = 0.3, box.padding = 0.25
  ) +
  
  # colorblind-friendly palette
  scale_color_manual(values = c(NS = "#9E9E9E", Down = "#D55E00", Up = "#0072B2")) +
  scale_alpha_identity() +
  
  # a touch of headroom above the tallest point
  coord_cartesian(ylim = c(0, y_max * 1.05)) +
  
  # crisp print style; no axis titles (you'll add them later)
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks   = element_line(color = "black"),
    axis.line    = element_line(color = "black"),
    axis.title   = element_blank()
  )

# Preview
print(p_volcano)

# 12) Save outputs (vector PDF + high-res PNG) ---------------------------------
ggsave("volcano_IL8.pdf", p_volcano, width = 7, height = 5, device = cairo_pdf)
ggsave("volcano_IL8.png", p_volcano, width = 7, height = 5, dpi = 600)
