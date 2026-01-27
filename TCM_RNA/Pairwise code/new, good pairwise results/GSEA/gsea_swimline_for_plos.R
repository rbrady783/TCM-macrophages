# ==========================================
# PLOS ONE GSEA DOTPLOT - SWIMLANE LAYOUT
# ==========================================

library(ggplot2)
library(dplyr)
library(tidyr)

dat <- read.csv("select_gsea_results_for_figure.csv", stringsAsFactors = FALSE)

# Pathway name formatting with abbreviations
format_pathway_names <- function(pathway_names) {
  cleaned <- gsub("^(HALLMARK_|GSE\\d+_|MEK_UP\\.V1_|FOROUTAN_|REACTOME_|MISHRA_|WP_|PDGF_UP\\.V1_)", "", pathway_names)
  cleaned <- gsub("(_|\\b)(UP|DOWN|DN)(\\b|_)", " ", cleaned, ignore.case = TRUE)
  cleaned <- gsub("_+", " ", cleaned)
  cleaned <- gsub("\\s+", " ", cleaned)
  cleaned <- trimws(cleaned)
  
  is_empty <- nchar(cleaned) == 0 | cleaned %in% c("Up", "Dn", "up", "dn", "DOWN", "DOWN ")
  cleaned[is_empty] <- pathway_names[is_empty]
  
  t <- tools::toTitleCase(tolower(cleaned))
  t <- gsub("\\bDna\\b", "DNA", t)
  t <- gsub("\\bRna\\b", "RNA", t)
  t <- gsub("\\bTnf\\b", "TNF", t)
  t <- gsub("\\bIl\\b", "IL", t)
  t <- gsub("\\bVs\\b", "vs", t)
  
  # Abbreviations
  t <- gsub("Epithelial Mesenchymal Transition", "EMT", t)
  t <- gsub("Extracellular Matrix", "ECM", t)
  t <- gsub("Macrophage", "Mφ", t)
  t <- gsub("Monocyte", "Mono", t)
  t <- gsub("Organization", "Org", t)
  t <- gsub("Associated", "Assoc", t)
  t <- gsub("Fibroblast", "Fibro", t)
  t <- gsub("Degradation", "Degrad", t)
  t <- gsub("Proinflammatory", "Proinflamm", t)
  t <- gsub("Intratumoral", "Intratumor", t)
  
  return(t)
}

# Cytokine normalization
normalize_cytokine <- function(x) {
  case_when(
    tolower(x) == "tnf.a" | tolower(x) == "tnf-a" ~ "TNF-α",
    tolower(x) == "tgf.b" | tolower(x) == "tgf-b" ~ "TGF-β", 
    tolower(x) == "kc.like" | tolower(x) == "kc-like" ~ "KC-like",
    tolower(x) == "il.10" | tolower(x) == "il-10" ~ "IL-10",
    tolower(x) == "il.8" | tolower(x) == "il-8" ~ "IL-8",
    tolower(x) == "vegf" ~ "VEGF",
    tolower(x) == "ccl2" ~ "CCL2",
    TRUE ~ as.character(x)
  )
}

# Prepare data - FIXED gene ratio calculation
df <- dat %>%
  mutate(
    pathway_clean = format_pathway_names(Description),
    cytokine_clean = normalize_cytokine(Cytokine),
    qvalue = as.numeric(qvalue),
    
    # Count core enrichment genes from the list (separated by "/")
    core_count = sapply(strsplit(as.character(core_enrichment), "/"), length),
    
    # Calculate gene ratio correctly: core genes / total pathway size
    gene_ratio = core_count / setSize
  ) %>%
  filter(is.finite(qvalue), qvalue <= 0.25)

# Order pathways alphabetically
pathway_order_levels <- sort(unique(df$pathway_clean))

# Custom cytokine order
cytokine_order_levels <- c("VEGF", "KC-like", "IL-8", "TGF-β", "CCL2", "TNF-α", "IL-10")

# Apply ordering
df <- df %>%
  mutate(
    pathway_order = factor(pathway_clean, levels = pathway_order_levels),
    cytokine_order = factor(cytokine_clean, levels = cytokine_order_levels)
  )

# SWIMLANE PLOT
swimlane_plot <- ggplot(df, aes(x = pathway_order, y = cytokine_order)) +
  geom_point(aes(size = gene_ratio, color = qvalue), alpha = 0.85) +  # Changed to gene_ratio and qvalue
  
  scale_color_gradient(
    low = "#0000FF",   
    high = "#FF0000", 
    trans = "log10",
    name = "FDR\nq-value",  # Changed label
    breaks = c(1e-04, 1e-03, 1e-02, 1e-01),
    labels = c("1e-04", "1e-03", "1e-02", "1e-01")
  ) +
  
  scale_size_continuous(name = "Gene\nRatio", range = c(2, 5.5)) +  # Changed label
  
  labs(x = NULL, y = NULL) +
  
  theme_minimal(base_size = 11, base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                               size = 9, color = "grey20", lineheight = 0.9,
                               margin = margin(t = 2)),
    axis.text.y = element_text(size = 11, color = "grey20", face = "bold",
                               margin = margin(r = 8)),
    
    panel.grid.major.x = element_line(color = "grey92", size = 0.3),
    panel.grid.major.y = element_line(color = "grey85", size = 0.4),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey40", size = 0.6, fill = NA),
    panel.background = element_rect(fill = "white"),
    
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = "grey40", size = 0.5),
    legend.title = element_text(size = 10, face = "bold", color = "black"),
    legend.text = element_text(size = 9, color = "grey20"),
    legend.key = element_rect(fill = "white", color = NA),
    
    plot.margin = margin(t = 60, r = 25, b = 25, l = 90),  # Reduced top margin
    plot.background = element_rect(fill = "white", color = NA)
  )

print(swimlane_plot)

# Save files
ggsave("gsea_swimlane.tiff", swimlane_plot, 
       width = 19.05, height = 14, units = "cm", 
       dpi = 600, compression = "lzw", bg = "white")
