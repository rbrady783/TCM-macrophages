# ══════════════════════════════════════════════════════════════════════════════
# VOLCANO PLOTS - PLOS ONE COMPLIANT
# Clean version with fixed line segments and optimized label positioning
# ══════════════════════════════════════════════════════════════════════════════

library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)

# ── SETTINGS ──────────────────────────────────────────────────────────────────
plot_width  <- 3.5
plot_height <- 2.1
plot_dpi    <- 600

padj_threshold <- 0.05059
lfc_threshold  <- 1
top_n_genes    <- 5

# ── VOLCANO PLOT FUNCTION ─────────────────────────────────────────────────────
create_volcano_plos <- function(res_df, cytokine_name, label_seed = 42, ycap = 8) {
  
  cat(paste("Processing", cytokine_name, "...\n"))
  
  # Prepare data
  volcano_df <- res_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      Hugo_Symbol = ifelse(is.na(Hugo_Symbol) | Hugo_Symbol == "" | Hugo_Symbol == "NA", 
                           gene, Hugo_Symbol),
      negLog10padj = -log10(padj),
      direction = case_when(
        padj < padj_threshold & log2FoldChange >= lfc_threshold ~ "Up",
        padj < padj_threshold & log2FoldChange <= -lfc_threshold ~ "Down",
        TRUE ~ "NS"
      ),
      above_cap = negLog10padj > ycap,
      y_plot = pmin(negLog10padj, ycap)
    ) %>%
    filter(!is.na(Hugo_Symbol), Hugo_Symbol != "", Hugo_Symbol != "NA")
  
  # Get top genes to label
  top_up <- volcano_df %>%
    filter(direction == "Up") %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    head(top_n_genes)
  
  top_down <- volcano_df %>%
    filter(direction == "Down") %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    head(top_n_genes)
  
  genes_to_label <- bind_rows(top_up, top_down)
  
  # Find outlier for asterisk
  outlier <- volcano_df %>%
    filter(above_cap) %>%
    slice_max(negLog10padj, n = 1)
  
  # Split data into layers
  volcano_df <- volcano_df %>%
    mutate(is_labeled = gene %in% genes_to_label$gene)
  
  ns_points <- volcano_df %>% filter(direction == "NS")
  sig_unlabeled <- volcano_df %>% filter(direction != "NS", !is_labeled)
  
  cat(paste("  Up:", nrow(top_up), "| Down:", nrow(top_down), "\n"))
  
  # Create plot
  p <- ggplot(volcano_df, aes(x = log2FoldChange, y = y_plot)) +
    
    # Layer 1: NS points
    geom_point(data = ns_points, color = "grey80", size = 1.5, alpha = 0.3) +
    
    # Layer 2: Significant unlabeled points (more transparent!)
    geom_point(data = sig_unlabeled, aes(color = direction), size = 1.5, alpha = 0.5) +
    
    # Threshold lines
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
               linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_hline(yintercept = -log10(padj_threshold),
               linetype = "dashed", color = "grey40", linewidth = 0.5) +
    
    # Cytokine label (lower right)
    annotate("text", x = Inf, y = -Inf, label = cytokine_name,
             hjust = 1.1, vjust = -0.5, size = 4.2, fontface = "bold",
             family = "Arial") +
    
    # Asterisk for capped values
    {if(nrow(outlier) > 0) {
      geom_text(data = outlier, aes(label = "*"), 
                vjust = -0.3, size = 4.2, color = "black", fontface = "bold",
                family = "Arial")
    }} +
    
    # Layer 3: Highlighted labeled points
    {if(nrow(genes_to_label) > 0) {
      geom_point(data = genes_to_label, aes(color = direction),
                 size = 3, stroke = 0.9, shape = 21, fill = "white")
    }} +
    
    # Layer 4: Gene labels with visible segments
    {if(nrow(genes_to_label) > 0) {
      geom_text_repel(data = genes_to_label, 
                      aes(label = Hugo_Symbol, color = direction),
                      size = 3, fontface = "bold", family = "Arial",
                      box.padding = 0.5,
                      point.padding = 0.4,
                      segment.color = "grey20",
                      segment.size = 0.5,
                      max.overlaps = Inf,
                      min.segment.length = 0,
                      force = 3,
                      seed = label_seed)
    }} +
    
    # Colors
    scale_color_manual(values = c(Down = "#3366CC", Up = "#DC3912")) +
    
    # Axis labels
    labs(x = expression(log[2]*" Fold Change"),
         y = expression(-log[10]*" (adj. p-value)")) +
    
    # Theme
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10, family = "Arial", face = "bold"),
      axis.text = element_text(size = 9, family = "Arial", color = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 12, r = 12, b = 12, l = 12),
      text = element_text(family = "Arial")
    )
  
  return(list(plot = p, data = volcano_df, labels = genes_to_label))
}

# ── PROCESS ALL CYTOKINES ─────────────────────────────────────────────────────
file_paths <- c(
  "vegf_quartile_results.csv",
  "il.8_quartile_results.csv",
  "il.10_quartile_results.csv",
  "kc.like_quartile_results.csv",
  "tnf.a_quartile_results.csv",
  "tgf.b_quartile_results.csv",
  "ccl2_quartile_results.csv"
)

cytokine_names <- c("VEGF", "IL-8", "IL-10", "KC-like", "TNF-α", "TGF-β", "CCL2")

# Custom seeds - optimized for each plot's specific issues
custom_seeds <- list(
  "IL-10" = 147,      # NR1H4 up
  "IL-8" = 1234,       # CCL7 up, HBP1 left
  "VEGF" = 999,       # GJA8 up, CD34 down, ARHGEF10L right
  "CCL2" = 321,       # Multiple overlaps
  "TGF-β" = 2024,     # Multiple overlaps
  "TNF-α" = 654,      # RNF14B up
  "KC-like" = 42      # Fine - no change
)

# Generate plots
for (i in seq_along(file_paths)) {
  tryCatch({
    res_df <- read_csv(file_paths[i], show_col_types = FALSE)
    
    seed_val <- ifelse(!is.null(custom_seeds[[cytokine_names[i]]]), 
                       custom_seeds[[cytokine_names[i]]], 42)
    
    result <- create_volcano_plos(res_df, cytokine_names[i], label_seed = seed_val)
    p <- result$plot
    
    safe_name <- gsub("[^A-Za-z0-9_-]", "", cytokine_names[i])
    
    # Save EPS
    ggsave(paste0(safe_name, "_volcano.eps"), p,
           width = plot_width, height = plot_height, 
           units = "in", device = cairo_ps, bg = "white")
    
    # Save TIFF
    ggsave(paste0(safe_name, "_volcano.tiff"), p,
           width = plot_width, height = plot_height,
           units = "in", dpi = plot_dpi, 
           compression = "lzw", bg = "white")
    
    cat(paste("✓", cytokine_names[i], "saved\n\n"))
    
  }, error = function(e) {
    cat(paste("✗ Error:", cytokine_names[i], "-", e$message, "\n\n"))
  })
}

# ══════════════════════════════════════════════════════════════════════════════
# IL-8 VOLCANO PLOT ONLY - With manual nudge for CCL7
# ══════════════════════════════════════════════════════════════════════════════

library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)

# ── SETTINGS ──────────────────────────────────────────────────────────────────
plot_width  <- 3.5
plot_height <- 2.1
plot_dpi    <- 600

padj_threshold <- 0.05059
lfc_threshold  <- 1
top_n_genes    <- 5

# Manual nudges for IL-8
manual_nudges <- list(
  "CCL7" = c(nudge_x = 0, nudge_y = 0.8)  # Adjust nudge_y higher if needed (try 1.0, 1.5, etc.)
)

# ── VOLCANO PLOT FUNCTION ─────────────────────────────────────────────────────
create_il8_volcano <- function(res_df, manual_nudges = NULL, ycap = 8) {
  
  cat("Processing IL-8...\n")
  
  # Prepare data
  volcano_df <- res_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      Hugo_Symbol = ifelse(is.na(Hugo_Symbol) | Hugo_Symbol == "" | Hugo_Symbol == "NA", 
                           gene, Hugo_Symbol),
      negLog10padj = -log10(padj),
      direction = case_when(
        padj < padj_threshold & log2FoldChange >= lfc_threshold ~ "Up",
        padj < padj_threshold & log2FoldChange <= -lfc_threshold ~ "Down",
        TRUE ~ "NS"
      ),
      above_cap = negLog10padj > ycap,
      y_plot = pmin(negLog10padj, ycap)
    ) %>%
    filter(!is.na(Hugo_Symbol), Hugo_Symbol != "", Hugo_Symbol != "NA")
  
  # Get top genes to label
  top_up <- volcano_df %>%
    filter(direction == "Up") %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    head(top_n_genes)
  
  top_down <- volcano_df %>%
    filter(direction == "Down") %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    head(top_n_genes)
  
  genes_to_label <- bind_rows(top_up, top_down)
  
  # Mark genes for manual adjustment
  genes_to_label$manual_adjust <- FALSE
  if(!is.null(manual_nudges)) {
    for(gene_name in names(manual_nudges)) {
      idx <- which(genes_to_label$Hugo_Symbol == gene_name)
      if(length(idx) > 0) {
        genes_to_label$manual_adjust[idx] <- TRUE
      }
    }
  }
  
  # Split into auto and manual groups
  genes_auto <- genes_to_label %>% filter(!manual_adjust)
  genes_manual <- genes_to_label %>% filter(manual_adjust)
  
  # Find outlier for asterisk
  outlier <- volcano_df %>%
    filter(above_cap) %>%
    slice_max(negLog10padj, n = 1)
  
  # Split data into layers
  volcano_df <- volcano_df %>%
    mutate(is_labeled = gene %in% genes_to_label$gene)
  
  ns_points <- volcano_df %>% filter(direction == "NS")
  sig_unlabeled <- volcano_df %>% filter(direction != "NS", !is_labeled)
  
  cat(paste("  Up:", nrow(top_up), "| Down:", nrow(top_down), "\n"))
  
  # Create plot
  p <- ggplot(volcano_df, aes(x = log2FoldChange, y = y_plot)) +
    
    # Layer 1: NS points
    geom_point(data = ns_points, color = "grey80", size = 1.5, alpha = 0.3) +
    
    # Layer 2: Significant unlabeled points
    geom_point(data = sig_unlabeled, aes(color = direction), size = 1.5, alpha = 0.5) +
    
    # Threshold lines
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
               linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_hline(yintercept = -log10(padj_threshold),
               linetype = "dashed", color = "grey40", linewidth = 0.5) +
    
    # Cytokine label (lower right)
    annotate("text", x = Inf, y = -Inf, label = "IL-8",
             hjust = 1.1, vjust = -0.5, size = 4.2, fontface = "bold",
             family = "Arial") +
    
    # Asterisk for capped values
    {if(nrow(outlier) > 0) {
      geom_text(data = outlier, aes(label = "*"), 
                vjust = -0.3, size = 4.2, color = "black", fontface = "bold",
                family = "Arial")
    }} +
    
    # Layer 3: Highlighted labeled points
    {if(nrow(genes_to_label) > 0) {
      geom_point(data = genes_to_label, aes(color = direction),
                 size = 3, stroke = 0.9, shape = 21, fill = "white")
    }} +
    
    # Layer 4: Auto-positioned gene labels
    {if(nrow(genes_auto) > 0) {
      geom_text_repel(data = genes_auto, 
                      aes(label = Hugo_Symbol, color = direction),
                      size = 3, fontface = "bold", family = "Arial",
                      box.padding = 0.5,
                      point.padding = 0.4,
                      segment.color = "grey20",
                      segment.size = 0.5,
                      max.overlaps = Inf,
                      min.segment.length = 0,
                      force = 3,
                      seed = 42)
    }} +
    
    # Layer 5: Manually nudged gene labels
    {if(nrow(genes_manual) > 0) {
      gene_name <- genes_manual$Hugo_Symbol[1]
      nudge_vals <- manual_nudges[[gene_name]]
      geom_text_repel(data = genes_manual, 
                      aes(label = Hugo_Symbol, color = direction),
                      size = 3, fontface = "bold", family = "Arial",
                      nudge_x = nudge_vals["nudge_x"],
                      nudge_y = nudge_vals["nudge_y"],
                      box.padding = 0.5,
                      point.padding = 0.4,
                      segment.color = "grey20",
                      segment.size = 0.5,
                      max.overlaps = Inf,
                      min.segment.length = 0)
    }} +
    
    # Colors
    scale_color_manual(values = c(Down = "#3366CC", Up = "#DC3912")) +
    
    # Axis labels
    labs(x = expression(log[2]*" Fold Change"),
         y = expression(-log[10]*" (adj. p-value)")) +
    
    # Theme
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10, family = "Arial", face = "bold"),
      axis.text = element_text(size = 9, family = "Arial", color = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 12, r = 12, b = 12, l = 12),
      text = element_text(family = "Arial")
    )
  
  return(p)
}

# ── GENERATE IL-8 PLOT ────────────────────────────────────────────────────────
res_df <- read_csv("il.8_quartile_results.csv", show_col_types = FALSE)

p <- create_il8_volcano(res_df, manual_nudges = manual_nudges)

# Save
ggsave("IL-8_volcano.eps", p,
       width = plot_width, height = plot_height, 
       units = "in", device = cairo_ps, bg = "white")

ggsave("IL-8_volcano.tiff", p,
       width = plot_width, height = plot_height,
       units = "in", dpi = plot_dpi, 
       compression = "lzw", bg = "white")

cat("\n✓ IL-8 volcano plot saved with manual nudge for CCL7\n")
cat("  If CCL7 needs to move more, increase nudge_y (try 1.0, 1.5, 2.0)\n")