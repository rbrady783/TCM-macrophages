# Volcano Plot Grid Creator from Existing PNG Files
# Purpose: Combine individual volcano plot PNG files into a 2x4 grid layout

# Load required libraries
library(cowplot)
library(grid)
library(magick) 
library(ggplot2)

# Read in your PNG volcano plots at high density for crisp text
il10_plot <- image_read("volcano_IL-10_top10.png", density = 1200)
tnfa_plot <- image_read("volcano_TNF-_top10.png", density = 1200)
il8_plot <- image_read("volcano_IL-8_top10.png", density = 1200)
kclike_plot <- image_read("volcano_KC-like_top10.png", density = 1200)
ccl2_plot <- image_read("volcano_CCL2_top10.png", density = 1200)
tgfb_plot <- image_read("volcano_TGF-_top10.png", density = 1200)
vegf_plot <- image_read("volcano_VEGF_top10.png", density = 1200)

# Convert to ggplot objects
il10_gg <- ggdraw() + draw_image(il10_plot)
tnfa_gg <- ggdraw() + draw_image(tnfa_plot)
il8_gg <- ggdraw() + draw_image(il8_plot)
kclike_gg <- ggdraw() + draw_image(kclike_plot)
ccl2_gg <- ggdraw() + draw_image(ccl2_plot)
tgfb_gg <- ggdraw() + draw_image(tgfb_plot)
vegf_gg <- ggdraw() + draw_image(vegf_plot)

# Create empty plot for the 8th position (like the legend in your original)
empty_gg <- ggdraw()

# Create 2x4 grid: 2 across, 4 down with reduced spacing
# Row 1: IL-10, TNF-α
# Row 2: IL-8, KC-like  
# Row 3: CCL2, TGF-β
# Row 4: VEGF, empty
final_plot <- plot_grid(
  il10_gg, tnfa_gg,
  il8_gg, kclike_gg,
  ccl2_gg, tgfb_gg,
  vegf_gg, empty_gg,
  ncol = 2, nrow = 4,
  align = "hv",
  axis = "tb"  # Align top-bottom to reduce vertical spacing
)

# Add the centered x-axis title at bottom
#final_with_title <- plot_grid(
 # final_plot,
  #ggdraw() + draw_label("log₂(Fold Change)", 
   #                     fontface = "bold", size = 20),
  #ncol = 1, nrow = 2,
  #rel_heights = c(20, 1)
#)

ggsave("combined_volcano_final.eps", final_plot, 
       width = 14, height = 20, units = "in", dpi = 300, 
       device = "eps", bg = "white")

# Save at high DPI for crisp text with adjusted dimensions to reduce white space
ggsave("combined_volcano_final.pdf", final_plot, 
       width = 14, height = 20, units = "in", dpi = 1200, bg = "white")
ggsave("combined_volcano_final.png", final_plot, 
       width = 14, height = 20, units = "in", dpi = 1200, bg = "white")

# Print completion message
cat("Combined volcano plot grid saved as:\n")
cat("  - combined_volcano_final.pdf\n")
cat("  - combined_volcano_final.png\n")

