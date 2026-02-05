library(ggbreak)
library(cowplot)

# Base plotting function (without labels)
plot_cytokine_base <- function(cyt) {
  df <- summary_df %>%
    filter(cytokine == cyt, treatment != "Ctrl") %>%
    mutate(
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
      panel.background    = element_rect(fill = "transparent", color = NA),
      plot.background     = element_rect(fill = "transparent", color = NA),
      legend.background   = element_rect(fill = "transparent", color = NA),
      legend.position     = "none",
      axis.title.x        = element_blank(), 
      axis.title.y        = element_blank(),
      axis.text.y         = element_text(face = "bold"),
      plot.title          = element_blank(),
      panel.grid.major.x  = element_line(color = "grey90"),
      panel.grid.minor.x  = element_blank()
    ) 
}

#IL-10 range check
il10_data_check <- summary_df %>% 
  filter(cytokine == "il.10", treatment != "Ctrl") %>%
  summarise(
    min_val = min(mean_z - se_z, na.rm = TRUE),
    max_val = max(mean_z + se_z, na.rm = TRUE)
  )
print(il10_data_check)

il10_values <- summary_df %>% 
  filter(cytokine == "il.10", treatment != "Ctrl") %>%
  arrange(mean_z) %>%
  select(treatment, mean_z, se_z) %>%
  mutate(
    lower = mean_z - se_z,
    upper = mean_z + se_z
  )
print(il10_values)

# IL-10 with two breaks for better spacing
# IL-10 with single break - let ggbreak handle labels automatically
p_il10_base <- plot_cytokine_base("il.10")
p_il10_broken <- p_il10_base +
  scale_y_break(
    breaks = c(8, 38),
    scales = 0.3,
    space = 0.05
    # Remove ticklabels parameter - let ggbreak auto-generate
  )

# Save without labels to avoid duplication (add labels manually later)
ggsave(
  filename = "IL-10_broken.png",
  plot = p_il10_broken,
  width = 8,
  height = 6,
  units = "in",
  dpi = 1200,
  bg = "transparent",
  type = "cairo-png"  # Better transparency handling
)

ggsave(
  filename = "IL-10_broken.pdf",
  plot = p_il10_broken,
  width = 8,
  height = 6,
  units = "in",
  bg = "transparent",
  device = cairo_pdf
)

# TNF-α range check
tnfa_specific <- summary_df %>% 
  filter(cytokine == "tnf.a", treatment %in% c("Nike", "OSA8")) %>%
  select(treatment, mean_z, se_z) %>%
  mutate(
    lower = mean_z - se_z,
    upper = mean_z + se_z
  )
print(tnfa_specific)

# TNF-α with breaks - avoid cutting through data points
p_tnfa_base <- plot_cytokine_base("tnf.a")
p_tnfa_broken <- p_tnfa_base +
  scale_y_break(
    breaks = c(6, 8.5),   # Break between low data (~0-6) and high data (12+)
    scales = 0.3,
    space = 0.05,
    ticklabels = c(0, 2, 4, 6, 8.5, 12, 16, 20, 24, 28, 32)
  )

# Save without labels to avoid duplication (add labels manually later)
ggsave(
  filename = "TNFa_broken.png",
  plot = p_tnfa_broken,
  width = 8,
  height = 6,
  units = "in",
  dpi = 1200,
  bg = "transparent",
  type = "cairo-png"  # Better transparency handling
)

ggsave(
  filename = "TNFa_broken.pdf",
  plot = p_tnfa_broken,
  width = 8,
  height = 6,
  units = "in",
  bg = "transparent",
  device = cairo_pdf
)
