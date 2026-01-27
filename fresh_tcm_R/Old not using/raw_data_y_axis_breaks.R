# ── 0) PACKAGES ────────────────────────────────────────────────────
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)    # for comma() in y-axis
library(ggbreak)   # for broken axes

# ── 1) READ YOUR CSV ───────────────────────────────────────────────
raw <- read_csv("raw_with_ctrls.csv", show_col_types = FALSE)

# ── 2) PIVOT TO LONG ───────────────────────────────────────────────
long <- raw %>%
  rename(cell_line = treatment) %>%              # rename “treatment” → “cell_line”
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

# ── 4) PLOTTING FUNCTION WITH FAINT VERTICAL GRID LINES ───────────
plot_raw_cytokine <- function(mycyt) {
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
      # vertical grid lines under each x-tick:
      panel.grid.major.x = element_line(color = "grey90"),
      panel.grid.minor.x = element_blank(),
      # no horizontal grid lines
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      
      axis.text.x     = element_text(angle = 45, hjust = 1),
      plot.title      = element_text(face = "bold", hjust = 0.5),
      legend.position = "right"
    ) +
    labs(
      title = toupper(mycyt),
      x     = NULL,
      y     = "pg/mL"
    )
}

# ── 5) PRODUCE & SAVE ALL SEVEN RAW PLOTS ──────────────────────────
cyts <- c("CCL2","IL.8","KC.LIKE","IL.10","TNF.A","TGF.B","VEGF") %>%
  tolower()

for(cyt in cyts) {
  p <- plot_raw_cytokine(cyt)
  ggsave(
    filename = paste0(toupper(gsub("\\.", "-", cyt)), "_raw.png"),
    plot     = p,
    width    = 8,
    height   = 6,
    units    = "in",
    dpi      = 600,
    bg       = "transparent"
  )
}

# ── 6) IL-10 BROKEN-AXIS PREVIEW & SAVE ────────────────────────────
p_il10       <- plot_raw_cytokine("il.10")

# Find the actual max to include in ticklabels (optional)
actual_max <- max(
  filter(long, cytokine=="il.10")$pg_ml,
  na.rm = TRUE
)

p_il10_broken <- p_il10 +
  scale_y_break(
    breaks     = list(c(2500, 5500), c(8000, 19000)),
    scales     = c(0.1, 0.4, 100),
    space      = c(0.02, 0.02),
    ticklabels = c(0, 2500, 5500, 8000, actual_max)
  ) +
  # remove the right‐hand y axis:
  theme(
    axis.text.y.right  = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right  = element_blank()
  )

print(p_il10_broken) 

ggsave(
  "IL-10_raw_broken.png",
  plot   = p_il10_broken,
  width  = 8, height = 6, units = "in",
  dpi    = 600,
  bg     = "transparent"
)

# ── 7) TNF-α BROKEN-AXIS PREVIEW & SAVE ───────────────────────────
p_tnfa       <- plot_raw_cytokine("tnf.a")

# compute max so we can put the upper break just above
#max_z <- long %>%
  #filter(cytokine=="tnf.a") %>%
  #pull(pg_ml) %>%
 # max(na.rm = TRUE)

p_tnfa_broken <- p_tnfa +
  scale_y_break(
    breaks     = list(c(400, 800)),
    scales     = c(0.3, 100),
    space      = c(0.02),
    ticklabels = c(0, 1000, 1300)
  ) +
  # remove the right‐hand y axis:
  theme(
    axis.text.y.right  = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right  = element_blank()
  )

print(p_tnfa_broken)
ggsave(
  "TNFa_raw_broken.png",
  plot   = p_tnfa_broken,
  width  = 8, height = 6, units = "in",
  dpi    = 600,
  bg     = "transparent"
)
