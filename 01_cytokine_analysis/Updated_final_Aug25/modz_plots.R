library(cowplot)
library(grid)
library(gridExtra)
library(magick) 
library(ggplot2)

# Read in your PNG plots at high density for crisp text
il10_plot <- image_read("IL-10_broken.png", density = 1200)
tnfa_plot <- image_read("TNFa_broken.png", density = 1200)
il8_plot <- image_read("IL-8_modz_plot.png", density = 1200)
kclike_plot <- image_read("KC-like_modz_plot.png", density = 1200)
ccl2_plot <- image_read("CCL2_modz_plot.png", density = 1200)
tgfb_plot <- image_read("TGF-beta_modz_plot.png", density = 1200)
vegf_plot <- image_read("VEGF_modz_plot.png", density = 1200)
legend_plot <- image_read("legend_histology.png", density = 1200)

# Convert to ggplot objects
il10_gg <- ggdraw() + draw_image(il10_plot)
tnfa_gg <- ggdraw() + draw_image(tnfa_plot)
il8_gg <- ggdraw() + draw_image(il8_plot)
kclike_gg <- ggdraw() + draw_image(kclike_plot)
ccl2_gg <- ggdraw() + draw_image(ccl2_plot)
tgfb_gg <- ggdraw() + draw_image(tgfb_plot)
vegf_gg <- ggdraw() + draw_image(vegf_plot)
# Make legend bigger by cropping/scaling
legend_gg <- ggdraw() + draw_image(legend_plot, scale = 1.3)

# Create 2x4 grid: 2 across, 4 down with reduced spacing
# Row 1: IL-10, TNF-α
# Row 2: IL-8, KC-like  
# Row 3: CCL2, TGF-β
# Row 4: VEGF, Legend
final_plot <- plot_grid(
  il10_gg, tnfa_gg,
  il8_gg, kclike_gg,
  ccl2_gg, tgfb_gg,
  vegf_gg, legend_gg,
  ncol = 2, nrow = 4,
  align = "hv",
  axis = "tb"  # Align top-bottom to reduce vertical spacing
)

# Add the centered x-axis title at bottom
final_with_title <- plot_grid(
  final_plot,
  ggdraw() + draw_label("Mean modified z-score (+/- SEM)", 
                        fontface = "bold", size = 20),
  ncol = 1, nrow = 2,
  rel_heights = c(20, 1)
)

# Save as EPS for journal submission
# 14 inches at 300 dpi = 4200 pixels (well above 789 minimum requirement)
# Aspect ratio maintained at 14:20 (0.7:1)
ggsave("combined_cytokine_final2.eps", final_with_title, 
       width = 14, height = 20, units = "in", dpi = 300, 
       device = "eps", bg = "white")

# Save at high DPI for crisp text with adjusted dimensions to reduce white space
ggsave("combined_cytokine_final2.pdf", final_with_title, 
       width = 14, height = 20, units = "in", dpi = 1200, bg = "white")

ggsave("combined_cytokine_final2.png", final_with_title, 
       width = 14, height = 20, units = "in", dpi = 1200, bg = "white")
