# ── LIBRARIES ─────────────────────────────────────────────────────────────────
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# ── 1) READ IN YOUR DATA ───────────────────────────────────────────────────────
treat_df <- read_csv("modz_with_histo.csv", show_col_types = FALSE)

# ── 2) SPECIFY YOUR CONTROL VALUES ─────────────────────────────────────────────
ctrl_df <- tibble::tribble(
  ~cytokine,  ~treatment, ~donor_1_modz,  ~donor_2_modz,  ~donor_3_modz,
  "vegf",     "Ctrl",      2.484984796,   -0.277050892,  -0.071416965,
  "il.8",     "Ctrl",     -1.908891365,   -0.518655989,   0.175225627,
  "kc.like",  "Ctrl",     -3.256047169,   -1.895588795,  -0.556340602,
  "il.10",    "Ctrl",     -0.775007357,   -0.479740681,  -0.619945737,
  "ccl2",     "Ctrl",     -1.110884361,    1.091099152,  -0.048639654,
  "tnf.a",    "Ctrl",     -0.673432753,   -0.682682226,  -0.674500000,
  "tgf.b",    "Ctrl",     -0.826084544,	  -0.685080499, 	0.078190664,
) %>%
  mutate(histology = "Control")

# ── 3) STACK & PIVOT TO LONG FORM ──────────────────────────────────────────────
full_df <- bind_rows(treat_df, ctrl_df)

long_df <- full_df %>%
  pivot_longer(
    cols         = starts_with("donor_"),
    names_to     = "donor",
    names_prefix = "donor_",
    values_to    = "z_mod"
  )

# ── 4) SUMMARIZE MEAN + SEM PER CYTOKINE × TREATMENT × HISTOLOGY ─────────────
summary_df <- long_df %>%
  group_by(cytokine, treatment, histology) %>%
  summarise(
    mean_z = mean(z_mod, na.rm = TRUE),
    se_z   = sd(z_mod,   na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# ── 5) EXTRACT CONTROL MEANS FOR THE DASHED LINE ───────────────────────────────
ctrl_line <- summary_df %>%
  filter(treatment == "Ctrl") %>%
  select(cytokine, ctrl_mean = mean_z)

# ── 6) DEFINE YOUR COLORS ─────────────────────────────────────────────────────
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

# ── 7) PLOTTING FUNCTION (POINTS + ERRORBARS + TRANSPARENT BACKGROUND) ────────
plot_cytokine <- function(cyt) {
  # map raw key to pretty title with Greek letters
  title_label <- switch(cyt,
                        "vegf"    = expression(bold(VEGF)),
                        "il.8"    = expression(bold(IL-8)),
                        "il.10"   = expression(bold(IL-10)),
                        "ccl2"    = expression(bold(CCL2)),
                        "kc.like" = expression(bold(KC~"-like")),
                        "tnf.a"   = expression(bold(TNF*alpha)),
                        "tgf.b"   = expression(bold(TGF*beta)),
                        cyt
  )
  
  df <- summary_df %>%
    filter(cytokine == cyt, treatment != "Ctrl") %>%
    mutate(
      # reorder so highest mean is at top when flipped
      treatment = fct_reorder(treatment, mean_z, .desc = FALSE)
    )
  
  y0 <- ctrl_line %>%
    filter(cytokine == cyt) %>%
    pull(ctrl_mean)
  
  ggplot(df, aes(x = treatment, y = mean_z, color = histology)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_z - se_z, ymax = mean_z + se_z),
                  width = 0.2, size = 0.5) +
    geom_hline(yintercept = y0, linetype = "dashed", color = "grey40") +
    scale_color_manual(values = hist_colors) +
    coord_flip() +
    theme_classic(base_size = 14, base_family = "sans") +
    theme(
      # make backgrounds transparent
      panel.background    = element_rect(fill = "transparent", color = NA),
      plot.background     = element_rect(fill = "transparent", color = NA),
      legend.background   = element_rect(fill = "transparent", color = NA),
      
      legend.position     = "none",
      axis.title.y        = element_blank(),
      axis.text.y         = element_text(face = "bold"),
      plot.title          = element_text(size = 16, hjust = 0.5),
      panel.grid.major.x  = element_line(color = "grey90"),
      panel.grid.minor.x  = element_blank()
    ) +
    labs(
      title = title_label,
      y     = "Mean modified z-score ± SEM"
    )
}

# ── 8) SAVE EACH PLOT AS A TRANSPARENT PNG AT 600 DPI ─────────────────────────
cytokines <- c("vegf","il.8","kc.like","il.10","ccl2","tnf.a","tgf.b")

for(cyt in cytokines) {
  p <- plot_cytokine(cyt)
  
  pretty <- switch(cyt,
                   "vegf"    = "VEGF",
                   "il.8"    = "IL-8",
                   #"il.10"   = "IL-10",
                   "ccl2"    = "CCL2",
                   "kc.like" = "KC-like",
                   #"tnf.a"   = "TNF-alpha",
                   "tgf.b"   = "TGF-beta"
  )
  
  ggsave(
    filename = paste0(pretty, "_modz_plot.png"),
    plot     = p,
    dpi      = 600,
    width    = 7,      # inches
    height   = 6,      # inches
    units    = "in",
    bg       = "transparent"
  )
}

