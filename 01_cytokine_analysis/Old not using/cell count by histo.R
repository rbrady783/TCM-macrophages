# 1. Load libraries
library(tidyverse)
library(car)      # for Levene’s test

# 2. Read your data
#    Assumes "cell_counts.csv" has columns: Cell_Line,Histology,Cell_Count
df <- read_csv("cell_counts.csv") %>%
  mutate(histology = factor(histology))

# 3. Summarize by histology
df %>% 
  group_by(histology) %>% 
  summarise(
    n    = n(),
    mean = mean(cell_count, na.rm = TRUE),
    sd   = sd(cell_count,   na.rm = TRUE)
  )

# 4. Test homogeneity of variance (ANOVA assumption)
leveneTest(cell_count ~ histology, data = df)

# 5. One‐way ANOVA
aov_res <- aov(cell_count ~ histology, data = df)
summary(aov_res)

# 5a. Tukey’s HSD post‐hoc (if ANOVA is significant)
TukeyHSD(aov_res)

# 6. Nonparametric alternative: Kruskal–Wallis
kruskal.test(cell_count ~ histology, data = df)

library(tidyverse)
library(ggpubr)

# compute max for annotation
y_max <- max(df$cell_count, na.rm = TRUE) * 1.05

# grab p-value separately (optional)
kw_p <- kruskal.test(cell_count ~ histology, data = df)$p.value

ggplot(df, aes(histology, cell_count, fill = histology)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.8, color = "black") +
  geom_jitter(width = 0.2, size = 2, alpha = 0.9, show.legend = FALSE, color = "black") +
  # put p-value at top center, using normalized coords:
  stat_compare_means(
    method       = "kruskal.test",
    label        = "p.format",
    label.x.npc  = "middle",   # horizontal center
    label.y.npc  = "top",      # vertical top
    vjust        = -0.5        # nudge it just above the panel
  ) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0.02, 0.15))  # extra 15% headroom
  ) +
  # greyscale fill instead of color brewer
  scale_fill_grey(start = 0.9, end = 0.5) +
  # ensure the box outlines stay black
  scale_color_manual(values = "black") +
  labs(
    title = "Cell Count by Tumor Histology",
    x     = "Histology",
    y     = "Cell Count"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position   = "none",
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(face = "bold", hjust = 0.5),
    axis.title        = element_text(color = "black"),
    axis.text         = element_text(color = "black")
  )
