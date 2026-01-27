# =========================================
# Publication-Ready Spearman Heatmap:
# Colored for all pairs, labels only on significant ones
# =========================================

# 1) Load required packages ----
# install.packages(c("Hmisc","reshape2","dplyr","stringr","ggplot2"))
library(Hmisc)      # rcorr()
library(reshape2)   # melt()
library(dplyr)      # data wrangling
library(stringr)    # str_replace()
library(ggplot2)    # plotting

# 2) Read & clean the data ----
data_raw <- read.csv(
  "heatmap.csv",
  row.names    = 1,
  check.names  = FALSE,
  fileEncoding = "UTF-8"
)

# 3) Tidy up column names ----
colnames(data_raw) <- colnames(data_raw) %>%
  str_replace("TNF-a", "TNF-α") %>%
  str_replace("TGF-b", "TGF-β") %>%
  gsub("\\.", "-", .)

# 4) Define display order of cytokines ----
analyte_order <- c("VEGF","IL-8","KC-like","IL-10","CCL2","TNF-α","TGF-β")
df_data       <- data_raw[, analyte_order]

# 5) Compute Spearman correlations & p-values ----
rc       <- rcorr(as.matrix(df_data), type = "spearman")
cor_mat  <- rc$r
p_mat    <- rc$P

# 6) Melt into long form & prepare two-line labels ONLY for p < 0.05 ----
merged_df <- left_join(
  melt(cor_mat,
       varnames   = c("Analyte1","Analyte2"),
       value.name = "rho"),
  melt(p_mat,
       varnames   = c("Analyte1","Analyte2"),
       value.name = "p_value"),
  by = c("Analyte1","Analyte2")
) %>%
  filter(Analyte1 != Analyte2) %>%   # drop diagonal
  mutate(
    Analyte1 = factor(Analyte1, levels = analyte_order),
    Analyte2 = factor(Analyte2, levels = analyte_order),
    signif   = case_when(
      p_value < 0.0001 ~ "****"
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    rho_lbl  = sprintf("%.2f", rho),
    # two-line label only if significant
    label2   = ifelse(p_value < 0.05,
                      paste0(rho_lbl, "\n", signif),
                      NA)
  )

sig_df <- merged_df %>% filter(!is.na(label2))

# 7) Plot all tiles, but text only on significant tiles ----
p <- ggplot() +
  # 7a) colored tiles for every pair
  geom_tile(
    data = merged_df,
    aes(x = Analyte1, y = Analyte2, fill = rho),
    color = "grey80", size = 0.2
  ) +
  # 7b) two-line labels only for significant pairs
  geom_text(
    data       = sig_df,
    aes(x = Analyte1, y = Analyte2, label = label2),
    color      = "black",
    size       = 4,
    lineheight = 0.8
  ) +
  # 7c) diverging color scale
  scale_fill_gradient2(
    low      = "#2166ac",
    mid      = "white",
    high     = "#b2182b",
    midpoint = 0,
    limits   = c(-1, 1),
    name     = expression(Spearman~rho)
  ) +
  # 7d) titles & theme
  labs(
    title    = "Pairwise Spearman Correlation of Cytokines",
    subtitle = "Only significant (p < 0.05) correlations are labeled"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title      = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle   = element_text(size = 12, hjust = 0.5, margin = margin(b = 8)),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y     = element_text(size = 12, color = "black"),
    legend.position = "right",
    legend.title    = element_text(face = "bold", size = 12),
    legend.text     = element_text(size = 10),
    panel.grid      = element_blank(),
    panel.border    = element_rect(color = "black", fill = NA, size = 0.5),
    axis.ticks      = element_line(color = "black")
  )

# 8) Display & save ----
print(p)
ggsave("correlation_heatmap_sig_labels_only.png",
       plot   = p,
       width  = 6,
       height = 5,
       dpi    = 600,
       bg     = "white")
