############################################################
# Cell counts by tumor histology — ANOVA 
# Modified for PLOS ONE figure requirements
############################################################

# ---- 0) Packages ----

library(tidyverse)  # readr/dplyr/ggplot2/tibble
library(car)        # Levene/Brown–Forsythe test
library(ggpubr)     # stat_compare_means for plotting p-values
library(scales)     # formatting axes

# load data

df <- read_csv(
  "cell_counts.csv",
  col_types = cols(
    cell_line   = col_character(),
    histo_5     = col_character(),
    histo_4     = col_character(),
    cell_count  = col_double()
  )
)

# 4-category histo
df4 <- df %>%
  transmute(
    cell_line,
    histology  = factor(histo_4),
    cell_count = cell_count
  )

# descriptive stats
sumtab <- df4 %>%
  group_by(histology) %>%
  summarise(
    n      = n(),
    mean   = mean(cell_count, na.rm = TRUE),
    sd     = sd(cell_count,   na.rm = TRUE),
    se     = sd / sqrt(n),
    median = median(cell_count, na.rm = TRUE),
    IQR    = IQR(cell_count,    na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(histology)

print(sumtab)

# Brown–Forsythe (Levene with center = median) for variance homogeneity
lev4 <- leveneTest(cell_count ~ histology, data = df4, center = median)
print(lev4)

# Fit one-way ANOVA
aov4 <- aov(cell_count ~ histology, data = df4)
print(summary(aov4))

# Diagnostics
plot(aov4, which = 1)
plot(aov4, which = 2)

# normality test of residuals
print(shapiro.test(residuals(aov4)))


# PLOS ONE compliant theme
# Using "sans" (system default sans-serif, similar to Arial)
# PLOS accepts Arial, Times, or Symbol - "sans" is universally available
# Using 11-12 pt for better legibility in publication
plos_theme <- theme_bw(base_size = 12, base_family = "sans") +
  theme(
    legend.position  = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.title       = element_text(color = "black", size = 11),
    axis.text        = element_text(color = "black", size = 10),
    plot.margin      = margin(2, 2, 2, 2, "pt")  # 2-point border
  )

# plot
p4 <- ggplot(df4, aes(histology, cell_count, fill = histology)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.8, color = "black") +
  geom_jitter(width = 0.2, size = 2, alpha = 0.9, show.legend = FALSE, 
              color = "black") +
  stat_compare_means(
    method      = "anova",
    label       = "p.format",
    label.x.npc = "middle",
    label.y.npc = "top",
    vjust       = -0.6,
    size        = 4.5  # Increased for better legibility (≈12 pt)
  ) +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0.02, 0.15))) +
  scale_fill_grey(start = 0.9, end = 0.5) +
  labs(
    x     = "Tumor Type",
    y     = "Cell Count (million)"
  ) +
  plos_theme

print(p4)

# 5-category encoding with OSA by itself
df5 <- df %>%
  transmute(
    cell_line,
    histology  = factor(histo_5),
    cell_count = cell_count
  )

# descriptive stats
sumtab_5 <- df5 %>%
  group_by(histology) %>%
  summarise(
    n      = n(),
    mean   = mean(cell_count, na.rm = TRUE),
    sd     = sd(cell_count,   na.rm = TRUE),
    se     = sd / sqrt(n),
    median = median(cell_count, na.rm = TRUE),
    IQR    = IQR(cell_count,    na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(histology)

print(sumtab_5)

# Brown–Forsythe (Levene with center = median) for variance homogeneity
lev5 <- leveneTest(cell_count ~ histology, data = df5, center = median)
print(lev5)

# Fit one-way ANOVA 
aov5 <- aov(cell_count ~ histology, data = df5)
print(summary(aov5))

# diagnostic plots
plot(aov5, which = 1)
plot(aov5, which = 2)

# normality test of residuals
print(shapiro.test(residuals(aov5)))


# plot
p5 <- ggplot(df5, aes(histology, cell_count, fill = histology)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.8, color = "black") +
  geom_jitter(width = 0.2, size = 2, alpha = 0.9, show.legend = FALSE, 
              color = "black") +
  stat_compare_means(
    method      = "anova",
    label       = "p.format",
    label.x.npc = "middle",
    label.y.npc = "top",
    vjust       = -0.8,
    size        = 4.5  # Increased for better legibility (≈12 pt)
  ) +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0.02, 0.15))) +
  scale_fill_grey(start = 0.9, end = 0.5) +
  labs(
    title = "Cell Count by Tumor Type",
    x     = "Tumor Type",
    y     = "Cell Count (million)"
  ) +
  plos_theme

print(p5)

# ---- Stack plots for PLOS ONE ----

library(cowplot)
library(grid)

# PLOS ONE dimensions: max height 8.75 in, width 7 in for full page
# Adjusted padding for shared labels
left_pad   <- 20
bottom_pad <- 20

p5_shared <- p5 +
  labs(x = NULL, y = NULL) +
  theme(plot.margin = margin(t = 10, r = 15, b = bottom_pad, l = left_pad))

p4_shared <- p4 +
  labs(x = NULL, y = NULL) +
  theme(plot.margin = margin(t = 10, r = 15, b = bottom_pad, l = left_pad))

# stack vertically
stacked <- plot_grid(
  p5_shared,
  p4_shared,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 1)
)

# shared labels (font family controlled by theme, not draw_label)
final_combo <- ggdraw(stacked) +
  draw_label("Cell Count (million)",
             x = 0.015, y = 0.5, angle = 90,
             hjust = 0.2, vjust = 0.5, fontface = "bold",
             size = 11) +
  draw_label("Tumor Type",
             x = 0.5, y = 0.015,
             hjust = 0.5, vjust = 0.2, fontface = "bold",
             size = 11)

print(final_combo)

# ---- Save in PLOS ONE compliant format ----
# TIFF with LZW compression, 300-600 dpi
# Max dimensions: 7.5 x 8.75 inches
# Using 7 x 8.5 inches to stay within limits

outdir <- "figs"
dir.create(outdir, showWarnings = FALSE)

# Save as TIFF with LZW compression (PLOS ONE required format)
ggsave(
  file.path(outdir, "Fig1.tif"),
  final_combo,
  width = 7,           # inches (within 2.63-7.5 range)
  height = 8.5,        # inches (within 8.75 max)
  units = "in",
  dpi = 600,           # 300-600 dpi required
  device = "tiff",
  compression = "lzw"  # LZW compression required
)

# Optional: Save individual panels if needed
ggsave(
  file.path(outdir, "Fig1A.tif"),
  p5,
  width = 7,
  height = 4,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw"
)

ggsave(
  file.path(outdir, "Fig1B.tif"),
  p4,
  width = 7,
  height = 4,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw"
)

# ---- Effect size (optional) ----
library(effectsize)

eta <- eta_squared(aov4, partial = FALSE)
omega <- omega_squared(aov4, partial = FALSE)

print(eta)
print(omega)

# ---- Important notes for PLOS ONE submission ----
cat("\n=== PLOS ONE Figure Checklist ===\n")
cat("✓ Format: TIFF with LZW compression\n")
cat("✓ Resolution: 600 dpi (300-600 required)\n")
cat("✓ Dimensions: 7 x 8.5 inches (within limits)\n")
cat("✓ Font: Sans-serif 10-12 pt (PLOS accepts Arial, Times, or Symbol)\n")
cat("✓ Color mode: RGB (default in R)\n")
cat("✓ File naming: Fig1.tif (update to match your figure number)\n")
cat("\nRemember to:\n")
cat("- Place figure captions in manuscript, NOT in figure file\n")
cat("- Cite as 'Fig 1' in text\n")
cat("- Keep file size < 10 MB\n")
cat("\nNote: Using 'sans' font (system default sans-serif, similar to Arial)\n")
cat("Alternative: Change base_family to 'Times' if you prefer Times New Roman\n")