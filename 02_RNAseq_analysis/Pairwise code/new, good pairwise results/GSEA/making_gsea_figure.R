# ==========================================
# CLEAN GSEA DOTPLOT - Alternative Cytokine Indicators
# ==========================================

library(ggplot2)
library(dplyr)

# ---- Read and prepare data ----
dat <- read.csv("select_gsea_results_for_figure.csv", stringsAsFactors = FALSE)

# Better pathway name cleaning - fix "Up"/"Dn" issues
format_pathway_names <- function(pathway_names) {
  cleaned <- gsub("^(HALLMARK_|GSE\\d+_|MEK_UP\\.V1_|FOROUTAN_|REACTOME_|MISHRA_|WP_|PDGF_UP\\.V1_)", "", pathway_names)
  # Remove trailing Up/Dn that creates empty names
  cleaned <- gsub("\\s+(Up|Dn)\\s*$", "", cleaned, ignore.case = TRUE)
  cleaned <- gsub("(_|\\b)(UP|DOWN|DN)(\\b|_)", " ", cleaned, ignore.case = TRUE)
  cleaned <- gsub("_+", " ", cleaned)
  cleaned <- gsub("\\s+", " ", cleaned)
  cleaned <- trimws(cleaned)
  
  # Filter out names that are too short or just "Up"/"Dn"
  cleaned[nchar(cleaned) < 3] <- paste("Pathway", seq_along(cleaned[nchar(cleaned) < 3]))
  
  t <- tools::toTitleCase(tolower(cleaned))
  t <- gsub("\\bDna\\b", "DNA", t)
  t <- gsub("\\bRna\\b", "RNA", t)
  t <- gsub("\\bTnf\\b", "TNF", t)
  t <- gsub("\\bIl\\b", "IL", t)
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
  # Order by gene ratio for diagonal line
  arrange(gene_ratio) %>%
  mutate(pathway_clean = factor(pathway_clean, levels = unique(pathway_clean)))

# CLEAN OPTION 1: Standard dotplot with manual right-side cytokine labels
# Create a mapping handling multiple cytokines per pathway
cytokine_mapping <- df %>%
  select(pathway_clean, cytokine_clean) %>%
  distinct() %>%
  group_by(pathway_clean) %>%
  summarise(cytokines = paste(cytokine_clean, collapse = ", "), .groups = "drop") %>%
  arrange(match(pathway_clean, levels(df$pathway_clean)))

clean_simple <- ggplot(df, aes(x = gene_ratio, y = pathway_clean)) +
  geom_point(aes(size = count, color = p.adjust)) +
  
  # Add cytokine labels on the right side - further out
  geom_text(data = cytokine_mapping,
            aes(x = max(df$gene_ratio) * 1.08, y = pathway_clean, label = cytokines),
            hjust = 0, vjust = 0.5, size = 2.8, color = "darkblue", fontface = "bold") +
  
  # Blue to red color scale
  scale_color_gradient(
    low = "#0000FF",   
    high = "#FF0000", 
    trans = "log10",
    name = "p.adjust",
    breaks = c(1e-04, 1e-03, 1e-02),
    labels = c("1e-04", "1e-03", "1e-02")
  ) +
  
  # EVEN SMALLER dot size range for option 3
  scale_size_continuous(name = "Count", range = c(1, 4)) +
  
  scale_x_continuous(name = "GeneRatio", 
                     # More space for cytokine labels
                     limits = c(0, max(df$gene_ratio) * 1.35),
                     expand = c(0.02, 0)) +
  scale_y_discrete(name = NULL) +
  
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 9),
    # MOVE LEGENDS TO BOTTOM
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = 15),
    # More right margin for cytokine labels
    plot.margin = margin(10, 120, 10, 10)
  ) +
  ggtitle("Clean Dotplot with Cytokine Labels")

print(clean_simple)

# OPTION 2: Facet by cytokine (multiple panels)
option2 <- ggplot(df, aes(x = gene_ratio, y = pathway_clean)) +
  geom_point(aes(size = count, color = p.adjust)) +
  
  # Blue to red color scale (matches reference)
  scale_color_gradient(
    low = "#0000FF",   
    high = "#FF0000", 
    trans = "log10",
    name = "p.adjust"
  ) +
  
  scale_size_continuous(name = "Count", range = c(3, 10)) +
  scale_x_continuous(name = "GeneRatio") +
  scale_y_discrete(name = NULL) +
  
  # Separate panel for each cytokine
  facet_wrap(~ cytokine_clean, ncol = 3, scales = "free_y") +
  
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 8),
    axis.text.y = element_text(size = 7),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  ) +
  ggtitle("Option 2: Separate panels by cytokine")

print(option2)

# OPTION 3: Add cytokine as suffix to pathway names
df_option3 <- df %>%
  mutate(pathway_with_cyto = paste0(pathway_clean, " (", cytokine_clean, ")")) %>%
  arrange(gene_ratio) %>%
  mutate(pathway_with_cyto = factor(pathway_with_cyto, levels = unique(pathway_with_cyto)))

option3 <- ggplot(df_option3, aes(x = gene_ratio, y = pathway_with_cyto)) +
  geom_point(aes(size = count, color = p.adjust)) +
  
  # Blue to red color scale
  scale_color_gradient(
    low = "#0000FF",   
    high = "#FF0000", 
    trans = "log10",
    name = "p.adjust"
  ) +
  
  scale_size_continuous(name = "Count", range = c(2, 6)) +
  scale_x_continuous(name = "GeneRatio") +
  scale_y_discrete(name = NULL) +
  
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 9),
    legend.position = "right"
  ) +
  ggtitle("")

print(option3)

# Save the clean options
ggsave("gsea_clean_standard.pdf", clean_simple, width = 10, height = 8)
ggsave("gsea_faceted_by_cytokine.pdf", option2, width = 14, height = 10)  
ggsave("gsea_cytokine_in_names.pdf", option3, width = 8, height = 10)

print("Three clean options created:")
print("1. Clean standard dotplot (no cytokine info shown)")
print("2. Faceted by cytokine (separate panels)") 
print("3. Cytokine name added to pathway labels")
print("Recommend option 2 (facets) or 3 (names in labels) for showing cytokine info")
