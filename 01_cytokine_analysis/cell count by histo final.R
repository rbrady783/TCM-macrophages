############################################################
# Cell counts by tumor histology — ANOVA 
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


#normality test of residuals
print(shapiro.test(residuals(aov4)))


# plot
p4 <- ggplot(df4, aes(histology, cell_count, fill = histology)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.8, color = "black") +
  geom_jitter(width = 0.2, size = 2, alpha = 0.9, show.legend = FALSE, 
              color = "black") +
  stat_compare_means(
    method      = "anova",     # show ANOVA p-value to align with the model
    label       = "p.format",
    label.x.npc = "middle",
    label.y.npc = "top",
    vjust       = -0.6
  ) +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0.02, 0.15))) +
  scale_fill_grey(start = 0.9, end = 0.5) +
  labs(
    x     = "Tumor Type",
    y     = "Cell Count (million)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position  = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    axis.title       = element_text(color = "black"),
    axis.text        = element_text(color = "black")
  )

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

#normality test of residuals
print(shapiro.test(residuals(aov5)))


# plot
p5 <- ggplot(df5, aes(histology, cell_count, fill = histology)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.8, color = "black") +
  geom_jitter(width = 0.2, size = 2, alpha = 0.9, show.legend = FALSE, 
              color = "black") +
  stat_compare_means(
    method      = "anova",     # show ANOVA p-value to align with the model
    label       = "p.format",
    label.x.npc = "middle",
    label.y.npc = "top",
    vjust       = -0.8
  ) +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0.02, 0.15))) +
  scale_fill_grey(start = 0.9, end = 0.5) +
  labs(
    title = "Cell Count by Tumor Type",
    x     = "Tumor Type",
    y     = "Cell Count (million)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position  = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    axis.title       = element_text(color = "black"),
    axis.text        = element_text(color = "black")
  )

print(p5)

# stack plots

library(cowplot)
library(grid)

# keep panels on the same y-scale 
#y_lim <- range(df$cell_count, na.rm = TRUE)

# extra padding (in points) so shared labels don't collide with ticks
left_pad   <- 20  # more space for shared Y label
bottom_pad <- 20  # more space for shared X label

p5_shared <- p5 +
  labs(x = NULL, y = NULL) +
 # coord_cartesian(ylim = y_lim) +
  theme(plot.margin = margin(t = 10, r = 15, b = bottom_pad, l = left_pad))

p4_shared <- p4 +
  labs(x = NULL, y = NULL) +
 # coord_cartesian(ylim = y_lim) +
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

# shared labels; nudge a bit farther from edges
final_combo <- ggdraw(stacked) +
  draw_label("Cell Count (million)",
             x = 0.015, y = 0.5, angle = 90,  # farther left
             hjust = 0.2, vjust = 0.5, fontface = "bold") +
  draw_label("Tumor Type",
             x = 0.5, y = 0.015,              # lower
             hjust = 0.5, vjust = 0.2, fontface = "bold")

print(final_combo)

# save (PNG @ 600 dpi + Cairo PDF)
outdir <- "figs"; dir.create(outdir, showWarnings = FALSE)
ggsave(file.path(outdir, "cell_counts_p5_p4_vertical.png"),
       final_combo, width = 7, height = 10, units = "in",
       dpi = 600, device = "png")
ggsave(file.path(outdir, "cell_counts_p5_p4_vertical.pdf"),
       final_combo, width = 7, height = 10, units = "in",
       device = grDevices::cairo_pdf)

###effect size- did not include
# Using effectsize package (easiest)
library(effectsize)

# Calculate both
eta <- eta_squared(aov4, partial = FALSE)
omega <- omega_squared(aov4, partial = FALSE)

print(eta)    # Biased estimate
print(omega)  # Unbiased estimate
