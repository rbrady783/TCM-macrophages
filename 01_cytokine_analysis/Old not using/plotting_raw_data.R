# ── 0) PACKAGES ───────────────────────────────────────────────────────────────
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)  # for comma() in y‐axis labels

# ── 1) READ YOUR CSV ──────────────────────────────────────────────────────────
raw <- read_csv("raw_with_ctrls.csv", show_col_types = FALSE)

# ── 2) PIVOT TO LONG ──────────────────────────────────────────────────────────
long <- raw %>%
  rename(cell_line = treatment) %>%                  
  pivot_longer(
    cols         = starts_with("donor_"),              
    names_to     = "donor",
    names_prefix = "donor_",
    values_to    = "pg_ml"
  ) %>%
  mutate(
    donor     = factor(donor, levels = c("1","2","3")),
    cell_line = factor(cell_line, levels = unique(cell_line))
  )

# ── 3) DEFINE YOUR DONOR COLORS ───────────────────────────────────────────────
donor_cols <- c(
  "1" = "#1f78b4",  # bright blue
  "2" = "#33a02c",  # bright green
  "3" = "#e31a1c"   # bright red
)

# ── 4) PLOTTING FUNCTION WITH VERTICAL GRID LINES ─────────────────────────────
plot_raw_cytokine <- function(mycyt) {
  df <- long %>% filter(cytokine == mycyt)
  
  ggplot(df, aes(x = cell_line, y = pg_ml, color = donor, group = donor)) +
    geom_line(size = 0.5, alpha = 0.3) +
    geom_point(size = 3) +
    scale_color_manual(
      values = donor_cols,
      labels = c("Dog 1","Dog 2","Dog 3"),
      name   = NULL
    ) +
    scale_y_continuous(labels = scales::comma) +
    theme_classic(base_size = 14) +
    theme(
      # turn on only the vertical grid under each x‐tick:
      panel.grid.major.x  = element_line(color = "grey90"),
      panel.grid.minor.x  = element_blank(),
      # remove horizontal grids
      panel.grid.major.y  = element_blank(),
      panel.grid.minor.y  = element_blank(),
      
      axis.text.x     = element_text(angle = 45, hjust = 1),
      strip.text      = element_text(face = "bold"),
      legend.position = "right",
      plot.title      = element_text(face = "bold", hjust = 0.5)
    ) +
    labs(
      title = toupper(mycyt),
      x     = NULL,
      y     = "pg/mL"
    )
}

# ── 5) LOOP OVER YOUR SEVEN CYTOKINES AND SAVE AT 600 DPI ───────────────────────
cyts <- c("vegf","il.8","kc.like","il.10","ccl2","tnf.a","tgf.b")

for(cyt in cyts) {
  p <- plot_raw_cytokine(cyt)
  ggsave(
    filename = paste0(toupper(gsub("\\.", "-", cyt)), "_raw.png"),
    plot     = p,
    width    = 8,
    height   = 5,
    units    = "in",
    dpi      = 600,
    bg       = "transparent"
  )
}
