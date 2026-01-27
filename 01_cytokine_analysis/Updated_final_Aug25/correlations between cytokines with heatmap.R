# =========================================
# Publication-Ready Spearman Heatmap:
# Colored for all pairs, labels only on significant ones
# WITH MULTIPLE COMPARISONS CORRECTION
# =========================================

# 1) Load required packages ----
# install.packages(c("Hmisc","reshape2","dplyr","stringr","ggplot2"))
#install.packages("Hmisc")
library(Hmisc)      # rcorr()
library(reshape2)   # melt()
library(dplyr)      # data wrangling
library(stringr)    # str_replace()
library(ggplot2)    # plotting

# 2) Read & clean the data ----
data_raw <- read.csv(
  "means_modz_for_heatmap.csv",
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

# *** NEW SECTION: Multiple Comparisons Correction *** ----
# Extract unique pairwise p-values (excluding diagonal and duplicates)
n_vars <- ncol(cor_mat)
upper_tri_indices <- which(upper.tri(p_mat), arr.ind = TRUE)
p_values_vector <- p_mat[upper_tri_indices]

# Apply FDR correction (Benjamini-Hochberg)
p_adj_vector <- p.adjust(p_values_vector, method = "fdr")

# Alternative: Bonferroni correction (more conservative)
#p_adj_vector <- p.adjust(p_values_vector, method = "bonferroni")

# Create adjusted p-value matrix
p_adj_mat <- matrix(1, nrow = n_vars, ncol = n_vars)
dimnames(p_adj_mat) <- dimnames(p_mat)

# Fill in the adjusted p-values (both upper and lower triangles)
p_adj_mat[upper_tri_indices] <- p_adj_vector
p_adj_mat[lower.tri(p_adj_mat)] <- t(p_adj_mat)[lower.tri(p_adj_mat)]

# Print detailed comparison of uncorrected vs corrected p-values
cat("Multiple Comparisons Summary:\n")
cat("Number of pairwise comparisons:", length(p_values_vector), "\n")
cat("Uncorrected significant (p < 0.05):", sum(p_values_vector < 0.05, na.rm = TRUE), "\n")
cat("FDR-corrected significant (p < 0.05):", sum(p_adj_vector < 0.05, na.rm = TRUE), "\n\n")

# Create a detailed comparison table for console viewing
comparison_df <- data.frame(
  Analyte1 = rep(rownames(cor_mat)[upper_tri_indices[,1]]),
  Analyte2 = rep(colnames(cor_mat)[upper_tri_indices[,2]]),
  Correlation = cor_mat[upper_tri_indices],
  P_uncorrected = p_values_vector,
  P_FDR_corrected = p_adj_vector,
  Sig_uncorrected = p_values_vector < 0.05,
  Sig_FDR = p_adj_vector < 0.05
) %>%
  arrange(P_uncorrected)

cat("Detailed P-value Comparison (sorted by uncorrected p-value):\n")
print(comparison_df, row.names = FALSE, digits = 4)
# *** END NEW SECTION *** ----

# 6) Melt into long form & prepare labels for BOTH corrected and uncorrected ----
merged_df <- left_join(
  melt(cor_mat,
       varnames   = c("Analyte1","Analyte2"),
       value.name = "rho"),
  melt(p_mat,
       varnames   = c("Analyte1","Analyte2"),
       value.name = "p_value"),
  by = c("Analyte1","Analyte2")
) %>%
  left_join(
    melt(p_adj_mat,
         varnames   = c("Analyte1","Analyte2"),
         value.name = "p_adjusted"),
    by = c("Analyte1","Analyte2")
  ) %>%
  filter(Analyte1 != Analyte2) %>%   # drop diagonal
  mutate(
    Analyte1 = factor(Analyte1, levels = analyte_order),
    Analyte2 = factor(Analyte2, levels = analyte_order),
    
    # Significance symbols for uncorrected p-values
    signif_raw = case_when(
      p_value < 0.0001 ~ "****",
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    
    # Significance symbols for FDR-corrected p-values
    signif_adj = case_when(
      p_adjusted < 0.0001 ~ "****",
      p_adjusted < 0.001 ~ "***",
      p_adjusted < 0.01  ~ "**",
      p_adjusted < 0.05  ~ "*",
      TRUE               ~ ""
    ),
    
    rho_lbl = sprintf("%.2f", rho),
    
    # Choose which correction to display (change this line to switch)
    # Option 1: Use FDR-corrected p-values (recommended for publication)
    label2 = ifelse(p_adjusted < 0.05,
                    paste0(rho_lbl, "\n", signif_adj),
                    NA),
    
    # Option 2: Use uncorrected p-values (for exploration)
    # label2 = ifelse(p_value < 0.05,
    #                 paste0(rho_lbl, "\n", signif_raw),
    #                 NA),
    
    # Option 3: Show both (for comparison)
    # label2 = ifelse(p_value < 0.05 | p_adjusted < 0.05,
    #                 paste0(rho_lbl, "\n", 
    #                        "Raw:", signif_raw, " FDR:", signif_adj),
    #                 NA)
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
    title    = " "
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

# 8) Display & save in multiple high-quality formats ----
print(p)

ggsave("correlation_heatmap_fdr_corrected.eps",
       plot   = p,
       width  = 6,
       height = 5,
       dpi    = 300,
       device = "eps",
       bg     = "white")

# Save as high-quality PNG
ggsave("correlation_heatmap_fdr_corrected.png",
       plot   = p,
       width  = 6,
       height = 5,
       dpi    = 600,
       bg     = "white")

# Save as high-quality PDF with Cairo
ggsave("correlation_heatmap_fdr_corrected.pdf",
       plot   = p,
       device = cairo_pdf,
       width  = 6,
       height = 5,
       bg     = "white")

# Create publication-ready supplementary table
supp_correlation_table <- comparison_df %>%
  mutate(
    Cytokine_1 = Analyte1,
    Cytokine_2 = Analyte2,
    `ρ` = sprintf("%.3f", Correlation),
    `p-value` = ifelse(P_uncorrected < 0.001, "<0.001", sprintf("%.3f", P_uncorrected)),
    `Adj. p-value` = ifelse(P_FDR_corrected < 0.001, "<0.001", sprintf("%.3f", P_FDR_corrected))
  ) %>%
  select(Cytokine_1, Cytokine_2, `ρ`, `p-value`, `Adj. p-value`) %>%
  arrange(Cytokine_1, Cytokine_2)

write.csv(supp_correlation_table, "Supplementary_Cytokine_Correlations.csv", row.names = FALSE)

# Print for verification
print(supp_correlation_table)
