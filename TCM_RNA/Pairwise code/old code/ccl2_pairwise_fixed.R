# =========================================
# ccl2 High vs. Low Pairwise DEG Analysis
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
ccl2_groups <- read.csv("ccl2_groups.csv", header = TRUE, 
                        stringsAsFactors = FALSE)

# 2) Subset your samples manually (as before) ----
raw_ccl2 <- raw %>%
  select(1,3,4,6,11,13,15,16,18:20,22,25)

# 3) Move gene IDs into rownames & drop that column ----
rownames(raw_ccl2) <- raw_ccl2[[1]]
raw_ccl2 <- raw_ccl2[, -1]

# 4) Check sample names against metadata ----
if (!all(colnames(raw_ccl2) %in% ccl2_groups$id)) {
  stop("Sample names in raw_ccl2 do not match ccl2_groups$id!")
} else {
  message("All sample names match.")
}

# 5) Prepare metadata ----
ccl2_groups$group <- factor(ccl2_groups$group, levels = c("low","high"))
rownames(ccl2_groups) <- ccl2_groups$id

# 6) Prefilter low‐count genes ----
count_threshold <- 3
min_samples <- 6  # e.g. half
keep <- rowSums(raw_ccl2 >= count_threshold) >= min_samples
filtered_raw_ccl2 <- raw_ccl2[keep, ]

# 7) Build & run DESeq2 ----
dds_ccl2 <- DESeqDataSetFromMatrix(
  countData = filtered_raw_ccl2,
  colData   = ccl2_groups[colnames(filtered_raw_ccl2), ],
  design    = ~ group
)

dds_ccl2 <- DESeq(dds_ccl2)

res_shrunk_ccl2 <- lfcShrink(
  dds_ccl2,
  coef = "group_high_vs_low",
  type = "ashr"      # <-- swap in ashr
)

# 9) Export results ----
res_df <- as.data.frame(res_shrunk_ccl2) %>%
  rownames_to_column("gene")
write.csv(res_df, "ccl2_pairwise_fixed.csv", row.names = FALSE)

# 9) Build a volcano data.frame -----------------------------------------------
# - Compute -log10(padj)
# - Classify each gene as Up / Down / NS at padj<0.05 and |log2FC|≥1
# - Make only the NS points semi-transparent to reduce "fog"
volcano_df2 <- res_df %>%
  mutate(
    negLog10padj = -log10(padj),
    direction    = case_when(
      padj < 0.05 & log2FoldChange >=  1 ~ "Up",
      padj < 0.05 & log2FoldChange <= -1 ~ "Down",
      TRUE                               ~ "NS"
    ),
    alpha_pt     = ifelse(direction == "NS", 0.25, 1)  # NS points are lighter
  ) %>%
  filter(!is.na(padj), !is.na(log2FoldChange))

# 10) Cap extreme -log10(padj) values so outliers don't dominate the y-scale ---
# - Everything above 'ycap' is drawn at ycap
# - We mark the single most extreme point with an asterisk on the plot
ycap <- 6
volcano_df2 <- volcano_df2 %>%
  mutate(
    above_cap = negLog10padj > ycap,
    y_plot    = ifelse(above_cap, ycap, negLog10padj)
  )

# pick the most extreme (if any) to annotate with an asterisk
outlier <- volcano_df2 %>%
  dplyr::filter(above_cap) %>%
  dplyr::slice_max(negLog10padj, n = 1)

# 11) (Optional) Highlight a gene of interest ----------------------------------
# - Replace 'ens.id' with any ENSCAFG you want to spotlight
# - Use the correct human-readable symbol (SNAI1, not "SNAIL1")
ens.id    <- "ENSCAFG00000011499"  # example: SNAI1
readable  <- "SNAI1"

highlight.df <- volcano_df2 %>%
  dplyr::filter(gene == ens.id) %>%
  mutate(symbol = readable)
# NOTE: If 'highlight.df' ends up empty (gene not in res), the highlight layers
#       will simply add nothing—no error.

# 12) Draw the volcano with high contrast + a "CCL2" panel tag -----------------
# - No axis titles (you'll add them later in layout/figure software)
# - Colorblind-friendly palette
# - Panel border and darker guide lines for print clarity
p_volcano <- ggplot(volcano_df2, aes(x = log2FoldChange, y = y_plot)) +
  # main scatter; only NS points semi-transparent
  geom_point(aes(color = direction, alpha = alpha_pt),
             size = 2.2, show.legend = FALSE) +
  
  # dashed significance thresholds (|log2FC|=1 and padj=0.05)
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed", color = "grey30", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "grey30", linewidth = 0.4) +
  
  # add a "CCL2" tag in the upper-left corner of the panel
  annotate("text", x = -Inf, y = Inf, label = "CCL2",
           hjust = -0.1, vjust = 1.2, size = 5, fontface = "bold") +
  
  # asterisk to indicate at least one value was capped at ycap (if present)
  geom_text(data = outlier, aes(label = "*"), vjust = -0.2, size = 5) +
  
  # highlight a chosen gene with a bordered point + label (if present)
  geom_point(
    data = highlight.df,
    shape = 21, fill = "gold", color = "black",
    size = 3.8, stroke = 0.6
  ) +
  ggrepel::geom_text_repel(
    data = highlight.df, aes(label = symbol),
    nudge_x = 0.6, nudge_y = 0.8, size = 4,
    segment.color = "black", segment.size = 0.3, box.padding = 0.25
  ) +
  
  # high-contrast, colorblind-friendly palette for directions
  scale_color_manual(values = c(NS = "#9E9E9E", Down = "red", Up = "green")) +
  scale_alpha_identity() +  # use the alpha values we computed (no legend)
  
  # keep a bit of headroom above the cap line
  coord_cartesian(ylim = c(0, ycap + 1)) +
  
  # classic theme + panel border for crisp print
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks   = element_line(color = "black"),
    axis.line    = element_line(color = "black"),
    axis.title   = element_blank()  # remove axis titles; you will add later
  )

# Preview in R
print(p_volcano)

# 13) Save outputs (vector PDF for manuscripts + high-res PNG) -----------------
# - 'cairo_pdf' produces very crisp text/lines in PDFs for journals
ggsave("volcano_ccl2_SNAI1.pdf", p_volcano, width = 7, height = 5, device = cairo_pdf)
ggsave("volcano_ccl2_SNAI1.png",  p_volcano, width = 7, height = 5, dpi = 600)
