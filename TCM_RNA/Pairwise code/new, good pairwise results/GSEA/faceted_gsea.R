# ==========================================
# PROFESSIONAL GSEA DOTPLOT - Faceted by Cytokine
# ==========================================

library(ggplot2)
library(dplyr)

# ---- Read and prepare data ----
dat <- read.csv("select_gsea_results_for_figure.csv", stringsAsFactors = FALSE)

# Improved pathway name cleaning with line breaks for long names
format_pathway_names <- function(pathway_names) {
  cleaned <- gsub("^(HALLMARK_|GSE\\d+_|MEK_UP\\.V1_|FOROUTAN_|REACTOME_|MISHRA_|WP_|PDGF_UP\\.V1_)", "", pathway_names)
  
  # Remove UP/DOWN/DN patterns but preserve the rest
  cleaned <- gsub("(_|\\b)(UP|DOWN|DN)(\\b|_)", " ", cleaned, ignore.case = TRUE)
  cleaned <- gsub("_+", " ", cleaned)
  cleaned <- gsub("\\s+", " ", cleaned)
  cleaned <- trimws(cleaned)
  
  # Only replace truly empty names
  is_empty <- nchar(cleaned) == 0 | cleaned %in% c("Up", "Dn", "up", "dn", "DOWN", "DOWN ")
  cleaned[is_empty] <- pathway_names[is_empty]
  
  t <- tools::toTitleCase(tolower(cleaned))
  t <- gsub("\\bDna\\b", "DNA", t)
  t <- gsub("\\bRna\\b", "RNA", t)
  t <- gsub("\\bTnf\\b", "TNF", t)
  t <- gsub("\\bIl\\b", "IL", t)
  t <- gsub("\\bVs\\b", "vs", t)
  
  # Add line breaks for long pathway names (over 40 characters)
  t <- sapply(t, function(name) {
    if (nchar(name) > 40) {
      # Find a good break point (space closest to middle)
      words <- strsplit(name, " ")[[1]]
      mid <- ceiling(length(words) / 2)
      line1 <- paste(words[1:mid], collapse = " ")
      line2 <- paste(words[(mid+1):length(words)], collapse = " ")
      return(paste(line1, line2, sep = "\n"))
    } else {
      return(name)
    }
  }, USE.NAMES = FALSE)
  
  return(t)
}

# Simple cytokine normalization
normalize_cytokine <- function(x) {
  case_when(
    x == "tnf.a" ~ "TNF-α",
    x == "tgf.b" ~ "TGF-β", 
    x == "kc.like" ~ "KC-like",
    x == "il.10" ~ "IL-10",
    x == "il.8" ~ "IL-8",
    x == "vegf" ~ "VEGF",
    x == "CCL2" ~ "CCL2",
    TRUE ~ toupper(x)
  )
}

# Prepare data
df <- dat %>%
  mutate(
    pathway_clean = format_pathway_names(Description),
    cytokine_clean = normalize_cytokine(Cytokine),
    p.adjust = as.numeric(p.adjust),
    count = setSize,
    gene_ratio = setSize / 500
  ) %>%
  filter(is.finite(p.adjust), p.adjust < 0.10) %>%
  arrange(gene_ratio) %>%
  mutate(pathway_clean = factor(pathway_clean, levels = unique(pathway_clean)))

# Professional faceted plot
professional_plot <- ggplot(df, aes(x = gene_ratio, y = pathway_clean)) +
  geom_point(aes(size = count, color = p.adjust), alpha = 0.8) +
  
  scale_color_gradient(
    low = "#0000FF",   
    high = "#FF0000", 
    trans = "log10",
    name = "Adjusted\np-value",
    breaks = c(1e-04, 1e-03, 1e-02),
    labels = c("1e-04", "1e-03", "1e-02")
  ) +
  
  scale_size_continuous(name = "Gene\nCount", range = c(2, 8)) +
  scale_x_continuous(name = "Gene Ratio") +
  scale_y_discrete(name = NULL) +
  
  facet_wrap(~ cytokine_clean, ncol = 3, scales = "free_y") +
  
  theme_bw(base_size = 13) +
  theme(
    # Darker, crisper text
    text = element_text(color = "black"),
    axis.text.x = element_text(color = "grey20", size = 12),
    axis.text.y = element_text(color = "grey20", size = 11, hjust = 1, lineheight = 0.9),
    axis.title = element_text(color = "black", size = 13, face = "bold"),
    
    # Cleaner grid
    panel.grid.major = element_line(color = "grey85", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey40", size = 0.8),
    
    # Professional facet labels
    strip.background = element_rect(fill = "grey90", color = "grey40", size = 0.8),
    strip.text = element_text(face = "bold", size = 13, color = "black"),
    
    # Legend styling
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = "grey40", size = 0.5),
    legend.title = element_text(size = 12, face = "bold", color = "black"),
    legend.text = element_text(size = 11, color = "grey20"),
    legend.key = element_rect(fill = "white", color = NA),
    
    # Margins - more space on left for pathway names
    plot.margin = margin(10, 15, 10, 10),
    plot.title = element_blank()
  )

print(professional_plot)

ggsave("gsea_professional_faceted.eps", professional_plot, 
       width = 16, height = 9, dpi = 300, device = cairo_ps, bg = "white")

# Save with better proportions - wider for 3 columns
ggsave("gsea_professional_faceted.pdf", professional_plot, width = 16, height = 9, device = cairo_pdf)
ggsave("gsea_professional_faceted.png", professional_plot, width = 16, height = 9, dpi = 600)

print("Professional faceted plot created!")
print("Improvements: Smaller dots (2-8 range), darker text, cleaner grid, professional styling")