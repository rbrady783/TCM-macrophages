library(dplyr)
library(stringr)
library(ggplot2)
library(readr)

# 1) read in your combined GSEA hits
#    assumes you have a data.frame `combined` with cols:
#      cytokine   (e.g. "IL6")
#      collection (one of "Hallmark","C2","C5","C7")
#      pathway    (the pathway name)
#      NES        (normalized enrichment score)
#      padj       (FDR adjusted p-value)
gsea_hits <- combined %>%
  filter(padj < 0.05)   # only keep significant

# 2) factor‐order your cytokines & pathways for nicer plotting
gsea_hits <- gsea_hits %>%
  mutate(
    cytokine = factor(cytokine, levels = unique(cytokine)),
    # order pathways by average NES
    pathway  = factor(pathway,
                      levels = gsea_hits %>%
                        group_by(pathway) %>%
                        summarize(avgNES = mean(NES)) %>%
                        arrange(avgNES) %>%
                        pull(pathway))
  )

# 3) bubble‐plot
p <- ggplot(gsea_hits, aes(x = cytokine, y = pathway)) +
  geom_point(aes(color = collection, size = abs(NES)), alpha = 0.75) +
  scale_color_brewer("Collection", palette = "Set2") +
  scale_size_continuous("|NES|", range = c(3, 8)) +
  labs(
    title = "GSEA hits across 7 cytokines × 4 collections",
    x     = "Cytokine",
    y     = "Pathway"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

print(p)

# 4) save it
ggsave("GSEA_bubbleplot.png", p, width = 8, height = 6, dpi = 300)
