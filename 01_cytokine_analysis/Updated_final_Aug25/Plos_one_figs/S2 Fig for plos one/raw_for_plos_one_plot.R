# ── 0) PACKAGES ────────────────────────────────────────────────────
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)    # for comma() in y-axis
library(ggbreak)   # for broken axes
library(cowplot)
library(grid)

# ── 1) READ YOUR CSV ───────────────────────────────────────────────
raw <- read_csv("raw_with_ctrls.csv", show_col_types = FALSE)

# ── 2) PIVOT TO LONG + CLEAN LABELS ───────────────────────────────
long <- raw %>%
  rename(cell_line = treatment) %>%
  # Clean whitespace from cell_line names to ensure uniform alignment
  mutate(cell_line = trimws(cell_line)) %>%
  pivot_longer(
    cols         = starts_with("donor_"),
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

# ── 4) BASE PLOTTING FUNCTION (PLOS ONE compliant) ────────────────
plot_raw_cytokine_base <- function(mycyt) {
  df <- filter(long, cytokine == mycyt)
  
  ggplot(df, aes(x = cell_line, y = pg_ml, color = donor, group = donor)) +
    geom_line(linewidth = 0.5, alpha = 0.3) +
    geom_point(size = 2) +
    scale_color_manual(
      values = donor_cols,
      labels = c("Dog 1","Dog 2","Dog 3"),
      name   = NULL
    ) +
    scale_y_continuous(
      labels = comma,
      expand = expansion(mult = c(0.02, 0.05))  # Padding for non-broken plots
    ) +
    theme_classic(base_size = 10, base_family = "sans") +
    theme(
      legend.position     = "none",
      
      # vertical grid lines under each x-tick:
      panel.grid.major.x = element_line(color = "grey90"),
      panel.grid.minor.x = element_blank(),
      # no horizontal grid lines
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      
      axis.text.x     = element_text(angle = 90, hjust = 1, vjust = 0.3, size = 7, color = "black"),
      axis.text.y     = element_text(size = 9, color = "black"),
      axis.title.x    = element_blank(),
      axis.title.y    = element_text(angle = 90, vjust = 0.5, size = 10, color = "black"),
      plot.title      = element_blank(),
      plot.margin     = unit(c(2, 2, 2, 2), "pt"),  # Use unit() instead of margin()
      text            = element_text(family = "sans")
    ) +
    labs(y = "pg/mL")
}

# ── 5) EXPORT SETTINGS (PLOS ONE compliant) ───────────────────────
# Wider plots for 4-across layout (gives more space for x-axis labels)
# ── 5) EXPORT SETTINGS (PLOS ONE compliant) ───────────────────────
# For 2 columns x 4 rows layout
plot_width  <- 3.4   # inches (2 across = 6.6", within 7.5" limit)
plot_height <- 2.1   # inches
plot_dpi    <- 600

# ── 6) SAVE REGULAR PLOTS ──────────────────────────────────────────
cytokines_regular <- c("vegf", "il.8", "kc.like", "ccl2", "tgf.b")

for (cyt in cytokines_regular) {
  p <- plot_raw_cytokine_base(cyt)
  
  # Create pretty name for file
  pretty <- switch(cyt,
                   "vegf"    = "VEGF",
                   "il.8"    = "IL-8",
                   "kc.like" = "KC-like",
                   "ccl2"    = "CCL2",
                   "tgf.b"   = "TGF-beta",
                   toupper(cyt))
  
  # Save as EPS (vector format)
  ggsave(paste0(pretty, "_raw.eps"), p,
         width = plot_width, height = plot_height,
         units = "in", device = cairo_ps, bg = "white")
  
  # Save as TIFF with LZW compression (PLOS requirement)
  ggsave(paste0(pretty, "_raw.tiff"), p,
         dpi = plot_dpi, width = plot_width, height = plot_height,
         units = "in", compression = "lzw", bg = "white")
}


# ── 9) CREATE AND SAVE LEGEND ─────────────────────────────────────
legend_data <- data.frame(
  donor = factor(c("1", "2", "3"), levels = c("1", "2", "3")),
  label = c("Dog 1", "Dog 2", "Dog 3")
)

# Create legend plot
legend_plot <- ggplot(legend_data, aes(x = 1, y = donor, color = donor)) +
  geom_point(size = 2) +
  scale_color_manual(
    values = donor_cols,
    labels = c("Dog 1", "Dog 2", "Dog 3"),
    name = "Donor"
  ) +
  guides(color = guide_legend(ncol = 1)) +
  theme_void() +
  theme(
    legend.title = element_text(size = 10, family = "sans", face = "bold"),
    legend.text  = element_text(size = 9, family = "sans"),
    legend.spacing.y = unit(0.15, "lines"),
    legend.key.size = unit(0.7, "lines"),
    text = element_text(family = "sans")
  )

# Extract legend
raw_legend <- get_legend(legend_plot)

# Save legend as TIFF (PLOS requirement)
tiff("legend_raw.tiff", width = plot_width, height = plot_height,
     units = "in", res = plot_dpi, compression = "lzw", bg = "white")
grid.newpage()
grid.draw(raw_legend)
dev.off()

# ── IL-10 STANDALONE ─────────────────────────────────────────────
il10_data <- long %>% filter(cytokine == "il.10")

p_il10_standalone <- ggplot(il10_data, aes(x = cell_line, y = pg_ml, color = donor, group = donor)) +
  geom_line(linewidth = 0.5, alpha = 0.3) +
  geom_point(size = 2) +
  scale_color_manual(
    values = donor_cols,
    labels = c("Dog 1","Dog 2","Dog 3"),
    name = NULL
  ) +
  scale_y_continuous(
    labels = comma,
    breaks = c(0, 500, 1000, 1500, 2000, 2500, 5500, 8000, 20000),  # Manual breaks
    limits = c(0, 20500),
    expand = c(0, 0)
  ) +
  scale_y_break(
    breaks = list(c(2500, 5500), c(8000, 19500)),
    scales = c(0.1, 0.1, 3),
    space = 0.2
  ) +
  theme_classic(base_size = 10, base_family = "sans") +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3, size = 7, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(angle = 90, vjust = 0.5, size = 10, color = "black"),
    plot.margin = unit(c(2, 2, 2, 2), "pt"),
    text = element_text(family = "sans")
  ) +
  labs(y = "pg/mL")

ggsave("IL10_raw_broken.tiff", p_il10_standalone,
       dpi = 600, width = 3.4, height = 2.1,
       units = "in", compression = "lzw", bg = "white")

# ── TNF-α STANDALONE ─────────────────────────────────────────────
tnfa_data <- long %>% filter(cytokine == "tnf.a")

p_tnfa_standalone <- ggplot(tnfa_data, aes(x = cell_line, y = pg_ml, color = donor, group = donor)) +
  geom_line(linewidth = 0.5, alpha = 0.3) +
  geom_point(size = 2) +
  scale_color_manual(
    values = donor_cols,
    labels = c("Dog 1","Dog 2","Dog 3"),
    name = NULL
  ) +
  scale_y_continuous(
    labels = comma,
    breaks = c(0, 100, 200, 300, 400, 800, 1200),  # Manual breaks
    limits = c(0, 1450),  # Hard-coded limits with headroom
    expand = c(0, 0)
  ) +
  scale_y_break(
    breaks = list(c(400, 800)),
    scales = c(0.2, 1),  # Keep bottom normal, compress top
    space = 0.2
  ) +
  theme_classic(base_size = 10, base_family = "sans") +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3, size = 7, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(angle = 90, vjust = 0.5, size = 10, color = "black"),
    plot.margin = unit(c(2, 2, 2, 2), "pt"),
    text = element_text(family = "sans")
  ) +
  labs(y = "pg/mL")

ggsave("TNFa_raw_broken.tiff", p_tnfa_standalone,
       dpi = 600, width = 3.4, height = 2.1,
       units = "in", compression = "lzw", bg = "white")
