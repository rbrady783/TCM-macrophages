library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

# 1) read your summary CSV
summary_df <- read_csv("summary_GSEA.csv", show_col_types = FALSE) %>%
  # ensure cytokine is a factor in the exact order you want
  mutate(
    cytokine = factor(cytokine, levels = c("CCL2", "KC-like", "TGFb", "TNFa")),
    # wrap only the long pathway names, not the cytokine labels
    pathway_short = str_wrap(pathway, width = 30)
  )

# 2) plot
p <- ggplot(summary_df, aes(x = cytokine, y = pathway_short)) +
  geom_point(aes(color = collection, size = NES), alpha = 0.9) +
  scale_x_discrete(labels = c(
    CCL2      = "CCL2",
    `KC-like` = "KC-like",
    TGFb      = expression(TGF~beta),
    TNFa      = expression(TNF~alpha)
  )) +
  scale_color_brewer("Collection", palette = "Set2") +
  scale_size_continuous(name = "NES", range = c(5,10)) +
  guides(
    color = guide_legend(override.aes = list(size = 6)),
    size  = guide_legend()
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x      = element_text(color = "black", size = 14,
                                    angle = 45, hjust = 1),
    axis.text.y      = element_text(color = "black", size = 14),
    axis.title       = element_blank(),
    legend.title     = element_text(size = 14, face = "bold"),
    legend.text      = element_text(size = 12),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) 

print(p)

ggsave(
  "GSEA_dotplot_selected_clean.png",
  p,
  width  = 8,
  height = 5,
  dpi    = 600
)
