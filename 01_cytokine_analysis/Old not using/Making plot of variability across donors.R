library(tidyverse)
library(lme4)      # if you also want the ICCs
library(RColorBrewer)

# assume df_long exists: cols = donor, cytokine, treatment, value

df_long <- df_long %>%
  mutate(cytokine = factor(cytokine)) %>%           # ensure it’s a factor
  mutate(cytokine = fct_recode(cytokine,
                               "VEGF"   = "vegf",
                               "IL-10"  = "il.10",
                               "IL-8"   = "il.8",
                               "TNF-α"  = "tnf.a",
                               "TGF-β"  = "tgf.b",
                               "KC-like"= "kc.like",
                               "CCL2"   = "ccl2"
  )) %>%
  # set the display order you want:
  mutate(cytokine = fct_relevel(cytokine,
                                "VEGF","IL-8","IL-10","TNF-α","TGF-β","KC-like","CCL2"
  ))


# 2. compute per‐cytokine mean ± SD across donors
cyto_stats <- df_long %>%
  group_by(cytokine, donor) %>%
  summarise(donor_mean = mean(value), .groups="drop") %>%
  group_by(cytokine) %>%
  summarise(
    mean_val = mean(donor_mean),
    sd_val   = sd(donor_mean),
    .groups  = "drop"
  )

# 3. choose a nice qualitative palette from RColorBrewer
palette_cols <- brewer.pal(n = 7, name = "Set2")

# 4. plot
ggplot(cyto_stats, aes(x = cytokine, y = mean_val, fill = cytokine)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                width = 0.2, size = 0.8) +
  scale_fill_manual(values = palette_cols) +
  labs(
    title = "Average Cytokine Expression Across Donors",
    subtitle = "Mean ± 1 SD (n = 3 donors)",
    x = NULL,
    y = "Normalized expression",
    caption = "Data normalized to modified z-score"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title     = element_text(face = "bold", size = 16),
    plot.subtitle  = element_text(size = 12, margin = margin(b = 10)),
    axis.text.x    = element_text(angle = 25, hjust = 1),
    axis.title.y   = element_text(margin = margin(r = 10)),
    legend.position = "none"
  )
