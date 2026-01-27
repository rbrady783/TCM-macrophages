# ──────────────────────────────────────────────────────────────────────────────
# ──────────────────────────────────────────────────────────────────────────────
# 0. If you don’t already have these:
# install.packages(c("readr","dplyr","tidyr","ggplot2","forcats"))
# ──────────────────────────────────────────────────────────────────────────────

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# 1) Read in treatment + histology data
treat_df <- read_csv("modz_with_histo.csv", show_col_types = FALSE)

# 2) Specify your control values
ctrl_df <- tibble::tribble(
  ~cytokine,  ~treatment, ~donor_1_modz,  ~donor_2_modz,  ~donor_3_modz,
  "vegf",     "Ctrl",      2.484984796,   -0.277050892,  -0.071416965,
  "il.8",     "Ctrl",     -1.908891365,   -0.518655989,   0.175225627,
  "kc.like",  "Ctrl",     -3.256047169,   -1.895588795,  -0.556340602,
  "il.10",    "Ctrl",     -0.775007357,   -0.479740681,  -0.619945737,
  "ccl2",     "Ctrl",     -1.110884361,    1.091099152,  -0.048639654,
  "tnf.a",    "Ctrl",     -0.673432753,   -0.682682226,  -0.674500000,
  "tgf.b",    "Ctrl",      0.974782498,    1.667203949,   2.300644251
) %>% 
  mutate(histology = "Control")

# 3) Stack & pivot to long form
full_df <- bind_rows(treat_df, ctrl_df)

long_df <- full_df %>%
  pivot_longer(
    cols         = starts_with("donor_"),
    names_to     = "donor",
    names_prefix = "donor_",
    values_to    = "z_mod"
  )

# 4) Summarize mean + SEM per cytokine × treatment × histology
summary_df <- long_df %>%
  group_by(cytokine, treatment, histology) %>%
  summarise(
    mean_z = mean(z_mod, na.rm = TRUE),
    se_z   = sd(z_mod,   na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# 5) Extract control means for dashed line
ctrl_line <- summary_df %>%
  filter(treatment == "Ctrl") %>%
  select(cytokine, ctrl_mean = mean_z)

# 6) Custom colors
hist_colors <- c(
  HS      = "indianred",
  HSA     = "navyblue",
  Lym     = "violetred",
  MC      = "lawngreen",
  Mel     = "black",
  OSA     = "dodgerblue",
  STS     = "mediumblue",
  TC      = "darkgreen",
  UC      = "seagreen1",
  Control = "#000000"
)

# 7) Plotting function using points + error bars, sans font, pretty titles
plot_cytokine <- function(cyt) {
  # A) map raw key to publication‐ready label (with Greek letters)
  title_label <- switch(cyt,
                        "vegf"    = "VEGF",
                        "il.8"    = "IL-8",
                        "il.10"   = "IL-10",
                        "ccl2"    = "CCL2",
                        "kc.like" = "KC-like",
                        "tnf.a"   = expression(TNF*alpha),
                        "tgf.b"   = expression(TGF*beta),
                        cyt       # fallback
  )
  
  # B) subset & reorder treatments by descending mean
  df <- summary_df %>%
    filter(cytokine == cyt, treatment != "Ctrl") %>%
    mutate(
      treatment = fct_reorder(treatment, mean_z, .desc = FALSE)
    )
  
  # C) control‐mean for this cytokine
  y0 <- ctrl_line %>%
    filter(cytokine == cyt) %>%
    pull(ctrl_mean)
  
  # D) build plot
  ggplot(df, aes(x = treatment, y = mean_z, color = histology)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_z - se_z, ymax = mean_z + se_z),
                  width = 0.2, size = 0.5) +
    geom_hline(yintercept = y0, linetype = "dashed", color = "grey40") +
    scale_color_manual(values = hist_colors) +
    coord_flip() +
    theme_classic(base_size = 14, base_family = "sans") +
    theme(
      legend.position     = "none",
      axis.title.y        = element_blank(),
      axis.text.y         = element_text(face = "bold"),
      plot.title          = element_text(face = "bold", size = 16, hjust = 0.5),
      panel.grid.major.x  = element_line(color = "grey90"),
      panel.grid.minor.x  = element_blank()
    ) +
    labs(
      title = title_label,
      y     = "Mean modified z-score ± SEM"
    )
}

# 8) Draw one plot at a time:
plot_cytokine("vegf")
plot_cytokine("il.8")
plot_cytokine("kc.like")
plot_cytokine("il.10")
plot_cytokine("ccl2")
plot_cytokine("tnf.a")
plot_cytokine("tgf.b")

# Vector of the cytokine codes you want to save
cytokines <- c("vegf","il.8","kc.like","il.10","ccl2","tnf.a","tgf.b")

# Loop over each, generate the plot, and save to disk
for(cyt in cytokines) {
  # 1) make the plot
  p <- plot_cytokine(cyt)
  
  # 2) Define a pretty filename
  pretty <- switch(cyt,
                   #"vegf"    = "VEGF",
                   #"il.8"    = "IL-8",
                   #"il.10"   = "IL-10",
                   #"ccl2"    = "CCL2",
                   #"kc.like" = "KC-like",
                   #"tnf.a"   = "TNF-alpha",
                   "tgf.b"   = "TGF-beta"
  )
  filename <- paste0(pretty, "_modz_plot.png")
  
  # 3) Save at 300 dpi, 6×4 inches
  ggsave(
    filename = filename,
    plot     = p,
    dpi      = 600,
    width    = 7,
    height   = 6,
    units    = "in",
    bg = "transparent"
  )
}

##Go to y_axis_breaks for adjusting of IL-10 and TNF-a