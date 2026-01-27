# 1) Make sure you’ve loaded ggbreak at the top of your script:
# install.packages("ggbreak")   # only once
library(ggbreak)

p2_il10 <- plot_cytokine("il.10")

p2_il10_broken <- p2_il10 +
  scale_y_break(
    breaks    = list(c(15, 23), c(25, 35), c(41,75)),
    scales    = c(0.2, 0.2, 0.2, 200),             # 3 segments: [0–15], [23–25], [38–80]
    space     = c(0.02, 0.02, 0.02),           # sliver for each break
    ticklabels = c(0, 5, 15, 25, 40)
  ) 

ggsave(
  filename = "IL-10_broken.png",
  plot     = p2_il10_broken,   # <- note: use the broken‐axis version
  width    = 7,
  height   = 6,
  units    = "in",
  dpi      = 600,
  bg       = "transparent" 
)

p2_tnfa <- plot_cytokine("tnf.a")

max_z <- max(summary_df %>% 
               filter(cytokine=="tnf.a") %>% pull(mean_z), na.rm=TRUE)

p2_tnfa_broken <- p2_tnfa +
  scale_y_break(
    list(c(8, 11), c(23, max_z + 1)),  # collapse 25→(just beyond your real max)
    scales    = c(0.2, 0.2, 1),        # now three panels, but 3rd is trivial
    space     = c(0.02, 0.02),
    ticklabels= c(0, 5, 8, 11, 25, round(max_z))
  )

print(p2_tnfa_broken)

ggsave(
  filename = "TNFa_broken.png",
  plot     = p2_tnfa_broken,   # <- note: use the broken‐axis version
  width    = 7,
  height   = 6,
  units    = "in",
  dpi      = 600,
  bg       = "transparent" 
)

##Legend
library(ggplot2)

# your colors
hist_colors <- c(
  OSA     = "dodgerblue",
  HSA     = "navyblue",
  STS     = "mediumblue",
  Mel     = "black",
  TC      = "darkgreen",
  MC      = "lawngreen",
  UC      = "seagreen1",
  HS      = "indianred",
  Lym     = "violetred"
)

hist_labels <- c(
  "Osteosarcoma",
  "Soft Tissue Sarcoma",
  "Hemangiosarcoma",
  "Melanoma",
  "Thyroid Carcinoma",
  "Mammary Carinoma",
  "Urothelial Carcinoma",
  "Histiocytic Sarcoma",
  "Lymphoma"
)

# dummy data: one point per histology level
legend_df <- data.frame(
  histology = factor(names(hist_colors), levels = names(hist_colors)),
  x = 1,
  y = seq_along(hist_colors)
)

# build a ggplot whose *only* purpose is to draw the legend
legend_plot <- ggplot(legend_df, aes(x, y, color = histology)) +
  geom_point(size = 4) +
  scale_color_manual(
    values = hist_colors,
    labels = hist_labels,
    name   = "Histology"
  ) +
  guides(color = guide_legend(ncol = 1)) +
  theme_void() +
  theme(
    legend.title    = element_text(face = "bold"),
    legend.text     = element_text(size = 10),
    legend.key.size = unit(0.6, "lines")
  )

# now pull out the legend grob
library(cowplot)
hist_legend <- get_legend(legend_plot)

# and display it
grid::grid.newpage()
grid::grid.draw(hist_legend)


