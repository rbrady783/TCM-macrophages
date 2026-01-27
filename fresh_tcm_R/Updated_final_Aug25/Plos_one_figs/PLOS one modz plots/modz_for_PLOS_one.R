# ── LIBRARIES ────────────────────────────────────────────────────────────────
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(cowplot)
library(grid)

# ── 1) READ DATA ─────────────────────────────────────────────────────────────
treat_df <- read_csv("modz_with_histo.csv", show_col_types = FALSE)

ctrl_df <- tibble::tribble(
  ~cytokine,  ~treatment, ~donor_1_modz,  ~donor_2_modz,  ~donor_3_modz,
  "vegf",     "Ctrl",      2.484984796,   -0.277050892,  -0.071416965,
  "il.8",     "Ctrl",     -1.908891365,   -0.518655989,   0.175225627,
  "kc.like",  "Ctrl",     -3.256047169,   -1.895588795,  -0.556340602,
  "il.10",    "Ctrl",     -0.775007357,   -0.479740681,  -0.619945737,
  "ccl2",     "Ctrl",     -1.110884361,    1.091099152,  -0.048639654,
  "tnf.a",    "Ctrl",     -0.673432753,   -0.682682226,  -0.6745,
  "tgf.b",    "Ctrl",     -0.826084544,   -0.685080499,   0.078190664
) %>% mutate(histology = "Control")

full_df <- bind_rows(treat_df, ctrl_df)

long_df <- full_df %>%
  pivot_longer(starts_with("donor_"), names_to = "donor", values_to = "z_mod")

summary_df <- long_df %>%
  group_by(cytokine, treatment, histology) %>%
  summarise(mean_z = mean(z_mod, na.rm = TRUE),
            se_z = sd(z_mod, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

ctrl_line <- summary_df %>%
  filter(treatment == "Ctrl") %>%
  select(cytokine, ctrl_mean = mean_z)

# ── 2) COLORS ────────────────────────────────────────────────────────────────
hist_colors <- c(
  HS = "indianred", HSA = "navyblue", Lym = "violetred",
  MC = "lawngreen", Mel = "black", OSA = "dodgerblue",
  STS = "mediumblue", TC = "darkgreen", UC = "seagreen1",
  Control = "#000000"
)

# ── 3) BASE PLOT ─────────────────────────────────────────────────────────────
plot_cytokine_base <- function(cyt) {
  df <- summary_df %>%
    filter(cytokine == cyt, treatment != "Ctrl") %>%
    mutate(treatment = fct_reorder(treatment, mean_z, .desc = FALSE))
  
  y0 <- ctrl_line |> filter(cytokine == cyt) |> pull(ctrl_mean)
  
  ggplot(df, aes(y = treatment, x = mean_z, color = histology)) +
    geom_point(size = 1.8) +
    geom_errorbarh(aes(xmin = mean_z - se_z, xmax = mean_z + se_z),
                   height = 0.25, linewidth = 0.5) +
    geom_vline(xintercept = y0, linetype = "dashed", color = "grey40") +
    scale_color_manual(values = hist_colors) +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 7, family = "Arial"),
      axis.text.x = element_text(size = 8, family = "Arial"),
      axis.ticks.length = unit(0.15, "cm"),
      axis.line.x = element_line(color = "black"),
      axis.line.x.top = element_blank(),
      axis.text.x.top = element_blank(),
      panel.grid.major.x = element_line(color = "grey90"),
      panel.grid.minor.x = element_blank(),
      plot.margin = margin(t = 8, r = 15, b = 10, l = 5),
      text = element_text(family = "Arial")
    )
}

# ── 4) EXPORT SETTINGS ──────────────────────────────────────────────────────
plot_width  <- 2.45
plot_height <- 2.85
plot_dpi    <- 600

cytokines_regular <- c("vegf", "il.8", "kc.like", "ccl2", "tgf.b")

for (cyt in cytokines_regular) {
  p <- plot_cytokine_base(cyt)
  pretty <- switch(cyt,
                   "vegf"    = "VEGF",
                   "il.8"    = "IL-8",
                   "kc.like" = "KC-like",
                   "ccl2"    = "CCL2",
                   "tgf.b"   = "TGF-beta",
                   cyt)
  
  ggsave(paste0(pretty, "_modz_plot.eps"), p,
         width = plot_width, height = plot_height,
         units = "in", device = cairo_ps, bg = "white")
  ggsave(paste0(pretty, "_modz_plot.tiff"), p,
         dpi = plot_dpi, width = plot_width, height = plot_height,
         units = "in", compression = "lzw", bg = "white")
}

# ── 5) LEGEND ───────────────────────────────────────────────────────────────
hist_data <- data.frame(
  code = c("OSA","STS","HSA","Mel","TC","MC","UC","HS","Lym"),
  label = c("Osteosarcoma","Soft Tissue Sarcoma","Hemangiosarcoma",
            "Melanoma","Thyroid Carcinoma","Mammary Carcinoma",
            "Urothelial Carcinoma","Histiocytic Sarcoma","Lymphoma/Leukemia")
)

hist_colors_legend <- setNames(
  c("dodgerblue","mediumblue","navyblue","black","darkgreen",
    "lawngreen","seagreen1","indianred","violetred"),
  hist_data$code
)

legend_df <- data.frame(histology = factor(hist_data$code, levels = hist_data$code),
                        x = 1, y = seq_along(hist_data$code))

p_legend <- ggplot(legend_df, aes(x, y, color = histology)) +
  geom_point(size = 1.8) +
  scale_color_manual(values = hist_colors_legend, labels = hist_data$label,
                     name = "Tumor Type") +
  guides(color = guide_legend(ncol = 1)) +
  theme_void() +
  theme(
    legend.title = element_text(size = 10, family = "Arial", face = "bold"),
    legend.text  = element_text(size = 9, family = "Arial"),
    legend.spacing.y = unit(0.15, "lines"),
    legend.key.size = unit(0.7, "lines"),
    text = element_text(family = "Arial")
  )

hist_legend <- get_legend(p_legend)

png("legend_histology.png", width = plot_width, height = plot_height,
    units = "in", res = plot_dpi, bg = "white")
grid.newpage(); grid.draw(hist_legend); dev.off()

cairo_ps("legend_histology.eps", width = plot_width, height = plot_height)
grid.newpage(); grid.draw(hist_legend); dev.off()

cat("\n✓ Regular plots + legend exported.\n")

################################################################
# ── LIBRARIES ────────────────────────────────────────────────────────────────
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggbreak)
library(cowplot)

# ── 1) READ DATA ─────────────────────────────────────────────────────────────
treat_df <- read_csv("modz_with_histo.csv", show_col_types = FALSE)

ctrl_df <- tibble::tribble(
  ~cytokine,  ~treatment, ~donor_1_modz,  ~donor_2_modz,  ~donor_3_modz,
  "vegf",     "Ctrl",      2.484984796,   -0.277050892,  -0.071416965,
  "il.8",     "Ctrl",     -1.908891365,   -0.518655989,   0.175225627,
  "kc.like",  "Ctrl",     -3.256047169,   -1.895588795,  -0.556340602,
  "il.10",    "Ctrl",     -0.775007357,   -0.479740681,  -0.619945737,
  "ccl2",     "Ctrl",     -1.110884361,    1.091099152,  -0.048639654,
  "tnf.a",    "Ctrl",     -0.673432753,   -0.682682226,  -0.6745,
  "tgf.b",    "Ctrl",     -0.826084544,   -0.685080499,   0.078190664
) %>% mutate(histology = "Control")

full_df <- bind_rows(treat_df, ctrl_df)

long_df <- full_df %>%
  pivot_longer(starts_with("donor_"), names_to = "donor", values_to = "z_mod")

summary_df <- long_df %>%
  group_by(cytokine, treatment, histology) %>%
  summarise(mean_z = mean(z_mod, na.rm = TRUE),
            se_z = sd(z_mod, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

ctrl_line <- summary_df %>%
  filter(treatment == "Ctrl") %>%
  select(cytokine, ctrl_mean = mean_z)

# ── 2) COLORS ────────────────────────────────────────────────────────────────
hist_colors <- c(
  HS = "indianred", HSA = "navyblue", Lym = "violetred",
  MC = "lawngreen", Mel = "black", OSA = "dodgerblue",
  STS = "mediumblue", TC = "darkgreen", UC = "seagreen1",
  Control = "#000000"
)

# ── 3) BROKEN PLOT FUNCTION ──────────────────────────────────────────────────
plot_broken_axis <- function(cyt, breaks_vec, ticks) {
  df <- summary_df %>%
    filter(cytokine == cyt, treatment != "Ctrl") %>%
    mutate(treatment = fct_reorder(treatment, mean_z, .desc = FALSE))
  
  y0 <- ctrl_line |> filter(cytokine == cyt) |> pull(ctrl_mean)
  
  ggplot(df, aes(y = treatment, x = mean_z, color = histology)) +
    geom_point(size = 1.8) +
    geom_errorbarh(aes(xmin = mean_z - se_z, xmax = mean_z + se_z),
                   height = 0.25, linewidth = 0.5) +
    geom_vline(xintercept = y0, linetype = "dashed", color = "grey40") +
    scale_color_manual(values = hist_colors) +
    scale_x_break(breaks = breaks_vec, scales = 0.3, space = 0.05) +
    scale_x_continuous(breaks = ticks) +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 7, family = "Arial"),
      axis.text.x = element_text(size = 8, family = "Arial"),
      axis.ticks.length = unit(0.15, "cm"),
      axis.line.x = element_line(color = "black"),
      axis.line.x.top = element_blank(),
      axis.text.x.top = element_blank(),
      axis.ticks.x.top = element_blank(),
      panel.grid.major.x = element_line(color = "grey92"),
      panel.grid.minor.x = element_blank(),
      plot.margin = margin(t = 8, r = 15, b = 10, l = 5),
      text = element_text(family = "Arial")
    )
}

# ── 4) EXPORT SETTINGS ──────────────────────────────────────────────────────
plot_width  <- 2.45
plot_height <- 2.85
plot_dpi    <- 600

# IL-10
p_il10 <- plot_broken_axis("il.10", breaks_vec = c(14, 30), ticks = c(0, 5, 10, 40, 60, 80))
ggsave("IL10_broken.eps", p_il10,
       width = plot_width, height = plot_height, units = "in", device = cairo_ps, bg = "white")
ggsave("IL10_broken.tiff", p_il10,
       dpi = plot_dpi, width = plot_width, height = plot_height,
       units = "in", compression = "lzw", bg = "white")

# TNF-a
p_tnfa <- plot_broken_axis("tnf.a", breaks_vec = c(6.2, 8.5), ticks = c(0, 2, 4, 6, 10, 20, 30))
ggsave("TNFa_broken.eps", p_tnfa,
       width = plot_width, height = plot_height, units = "in", device = cairo_ps, bg = "white")
ggsave("TNFa_broken.tiff", p_tnfa,
       dpi = plot_dpi, width = plot_width, height = plot_height,
       units = "in", compression = "lzw", bg = "white")

cat("\n✓ Broken-axis plots exported.\n")
