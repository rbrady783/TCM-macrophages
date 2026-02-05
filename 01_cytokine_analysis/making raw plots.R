# ── 0) PACKAGES ────────────────────────────────────────────────────
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)    # for comma() in y-axis
library(ggbreak)   # for broken axes
library(cowplot)
library(grid)
library(gridExtra)
library(magick)

# ── 1) READ YOUR CSV ───────────────────────────────────────────────
raw <- read_csv("raw_with_ctrls.csv", show_col_types = FALSE)

# ── 2) PIVOT TO LONG ───────────────────────────────────────────────
long <- raw %>%
  rename(cell_line = treatment) %>%              # rename "treatment" → "cell_line"
  pivot_longer(
    cols         = starts_with("donor_"),        # gather donor_1/2/3 only
    names_to     = "donor",
    names_prefix = "donor_",
    values_to    = "pg_ml"
  ) %>%
  mutate(
    donor     = factor(donor,     levels = c("1","2","3")),
    cell_line = factor(cell_line, levels = unique(cell_line))
  )

# ── 3) DEFINE YOUR DONOR COLORS ───────────────────────────────────
donor_cols <- c(
  "1" = "#1f78b4",  # bright blue
  "2" = "#33a02c",  # bright green
  "3" = "#e31a1c"   # bright red
)

# ── 4) BASE PLOTTING FUNCTION (without labels - you'll add manually) ───────────
plot_raw_cytokine_base <- function(mycyt) {
  df <- filter(long, cytokine == mycyt)
  
  ggplot(df, aes(x = cell_line, y = pg_ml, color = donor, group = donor)) +
    geom_line(size = 0.5, alpha = 0.3) +
    geom_point(size = 3) +
    scale_color_manual(
      values = donor_cols,
      labels = c("Dog 1","Dog 2","Dog 3"),
      name   = NULL
    ) +
    scale_y_continuous(labels = comma) +
    theme_classic(base_size = 14) +
    theme(
      panel.background    = element_rect(fill = "transparent", color = NA),
      plot.background     = element_rect(fill = "transparent", color = NA),
      legend.background   = element_rect(fill = "transparent", color = NA),
      legend.position     = "none",
      
      # vertical grid lines under each x-tick:
      panel.grid.major.x = element_line(color = "grey90"),
      panel.grid.minor.x = element_blank(),
      # no horizontal grid lines
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 12, color = "black"),  # Larger, darker
      axis.text.y     = element_text(size = 12, color = "black"),  # Larger, darker
      axis.title.x    = element_blank(),
      axis.title.y    = element_text(angle = 90, vjust = 0.5, size = 14, color = "black"),  # Larger, darker
      plot.title      = element_blank()
    ) +
    labs(y = "pg/mL")
}

# Regular plots WITH labels (for non-broken plots)
plot_raw_cytokine_with_label <- function(mycyt) {
  # Create proper cytokine label
  title_label <- switch(mycyt,
                        "vegf"    = "VEGF",
                        "il.8"    = "IL-8", 
                        "ccl2"    = "CCL2",
                        "kc.like" = "KC-like",
                        "tgf.b"   = "TGF-β",
                        toupper(mycyt)
  )
  
  plot_raw_cytokine_base(mycyt) +
    theme(
      axis.title.y = element_text(angle = 90, vjust = 0.5, size = 14, color = "black")
    ) +
    annotate("text", 
             x = Inf,
             y = Inf,
             label = title_label,
             hjust = 1.05,
             vjust = 1.2,
             size = 6,
             fontface = "bold")
}

# ── 5) PLOTTING FUNCTION WITH TITLES (for individual saves) ───────────
plot_raw_cytokine <- function(mycyt) {
  # Only add legend for VEGF plot, NO centered titles for any plots
  if(mycyt == "vegf") {
    plot_raw_cytokine_with_label(mycyt) +
      theme(
        axis.title.y = element_text(size = 12),
        legend.position = "right"  # Only show legend on VEGF
      )
  } else {
    plot_raw_cytokine_with_label(mycyt) +
      theme(
        axis.title.y = element_text(size = 12),
        legend.position = "none"  # No legend on other plots
      )
  }
}

# ── 6) SAVE INDIVIDUAL PLOTS AT HIGH DPI ──────────────────────────
cyts <- c("ccl2","il.8","kc.like","tgf.b","vegf")

for(cyt in cyts) {
  p <- plot_raw_cytokine_with_label(cyt)
  
  # Save PNG
  ggsave(
    filename = paste0(toupper(gsub("\\.", "-", cyt)), "_raw.png"),
    plot     = p,
    width    = 8,
    height   = 6,
    units    = "in",
    dpi      = 1200,
    bg       = "transparent"
  )
  
  # Save PDF with cairo
  ggsave(
    filename = paste0(toupper(gsub("\\.", "-", cyt)), "_raw.pdf"),
    plot     = p,
    width    = 8,
    height   = 6,
    units    = "in",
    bg       = "transparent",
    device   = cairo_pdf
  )
}

# ── 7) IL-10 WITH BROKEN AXIS (no label - you'll add manually) ──────────────
p_il10_base <- plot_raw_cytokine_base("il.10")

p_il10_broken <- p_il10_base +
  scale_y_break(
    breaks     = list(c(2500, 5500), c(8000, 19000)),
    scales     = c(0.1, 0.4, 100),
    space      = c(0.02, 0.02),
    ticklabels = c(0, 1000, 2000, 2500, 5500, 8000, 15000, 19000)  # Custom labels to avoid squishing
  ) +
  theme(
    axis.text.y.right  = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right  = element_blank()
  )

# Save broken axis plots with both PNG and PDF
ggsave("IL-10_raw_broken.png", p_il10_broken,
       width = 8, height = 6, units = "in", dpi = 1200, bg = "transparent")

ggsave("IL-10_raw_broken.pdf", p_il10_broken,
       width = 8, height = 6, units = "in", bg = "transparent", device = cairo_pdf)

# ── 8) TNF-α WITH BROKEN AXIS (no label - you'll add manually) ───────────────
p_tnfa_base <- plot_raw_cytokine_base("tnf.a")

p_tnfa_broken <- p_tnfa_base +
  scale_y_break(
    breaks     = list(c(400, 800)),
    scales     = c(0.3, 100),
    space      = c(0.02),
    ticklabels = c(0, 200, 400, 800, 1000, 1200)  # Custom labels to avoid squishing
  ) +
  theme(
    axis.text.y.right  = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right  = element_blank()
  )

ggsave("TNFa_raw_broken.png", p_tnfa_broken,
       width = 8, height = 6, units = "in", dpi = 1200, bg = "transparent")

ggsave("TNFa_raw_broken.pdf", p_tnfa_broken,
       width = 8, height = 6, units = "in", bg = "transparent", device = cairo_pdf)

# ── 9) CREATE LEGEND ────────────────────────────────────────────
# Create a dummy plot to extract legend
legend_plot <- ggplot(long, aes(x = cell_line, y = pg_ml, color = donor)) +
  geom_point(size = 4) +
  scale_color_manual(
    values = donor_cols,
    labels = c("Dog 1", "Dog 2", "Dog 3"),
    name = "Donor"
  ) +
  theme_void() +
  theme(
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    legend.spacing.y = unit(0.2, "lines"),
    legend.key.size = unit(0.8, "lines"),
    legend.background = element_rect(fill = "transparent", color = NA)
  )

# Extract and save legend
raw_legend <- get_legend(legend_plot)

png("legend_raw.png", width = 3, height = 4, units = "in", 
    res = 1200, bg = "transparent")
grid::grid.newpage()
grid::grid.draw(raw_legend)
dev.off()

# ── 10) COMBINE ALL PLOTS IN 2x4 LAYOUT ───────────────────────────
# Create separate legend for bottom-right position
legend_plot_combined <- ggplot(long, aes(x = cell_line, y = pg_ml, color = donor)) +
  geom_point(size = 4) +
  scale_color_manual(
    values = donor_cols,
    labels = c("Dog 1", "Dog 2", "Dog 3"),
    name = "Donor"
  ) +
  theme_void() +
  theme(
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    legend.spacing.y = unit(0.2, "lines"),
    legend.key.size = unit(0.8, "lines"),
    legend.background = element_rect(fill = "transparent", color = NA)
  )

# Extract legend and save
combined_legend <- get_legend(legend_plot_combined)

png("legend_raw_combined.png", width = 3, height = 4, units = "in", 
    res = 1200, bg = "transparent")
grid::grid.newpage()
grid::grid.draw(combined_legend)
dev.off()

# Read all images at high density
il10_plot <- image_read("IL-10_raw_broken.png", density = 1200)
tnfa_plot <- image_read("TNFa_raw_broken.png", density = 1200)
il8_plot <- image_read("IL-8_raw.png", density = 1200)
kclike_plot <- image_read("KC-LIKE_raw.png", density = 1200)
ccl2_plot <- image_read("CCL2_raw.png", density = 1200)
tgfb_plot <- image_read("TGF-B_raw.png", density = 1200)
vegf_plot <- image_read("VEGF_raw.png", density = 1200)  # Use regular VEGF without legend
legend_plot_img <- image_read("legend_raw_combined.png", density = 1200)

# Convert to ggplot objects
il10_gg <- ggdraw() + draw_image(il10_plot)
tnfa_gg <- ggdraw() + draw_image(tnfa_plot)
il8_gg <- ggdraw() + draw_image(il8_plot)
kclike_gg <- ggdraw() + draw_image(kclike_plot)
ccl2_gg <- ggdraw() + draw_image(ccl2_plot)
tgfb_gg <- ggdraw() + draw_image(tgfb_plot)
vegf_gg <- ggdraw() + draw_image(vegf_plot)
legend_gg <- ggdraw() + draw_image(legend_plot_img, scale = 1.3)

# Create 2x4 grid layout
final_plot <- plot_grid(
  il10_gg, tnfa_gg,
  il8_gg, kclike_gg,
  ccl2_gg, tgfb_gg,
  vegf_gg, legend_gg,
  ncol = 2, nrow = 4,
  align = "hv"
)

# ── 11) ADD FINAL LAYOUT WITHOUT BOTTOM TITLE ─────────────────────────
# Remove the bottom pg/mL title - just use the final plot as-is
final_plot_clean <- final_plot

# Save final combined plot without bottom title
ggsave("combined_raw_cytokine_final.pdf", final_plot_clean, 
       width = 14, height = 20, units = "in", dpi = 1200, bg = "white")

ggsave("combined_raw_cytokine_final.png", final_plot_clean, 
       width = 14, height = 20, units = "in", dpi = 1200, bg = "white")
