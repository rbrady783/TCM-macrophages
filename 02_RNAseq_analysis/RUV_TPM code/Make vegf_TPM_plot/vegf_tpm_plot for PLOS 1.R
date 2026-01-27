#— Install once if needed:
# install.packages(c("readr","tidyr","dplyr","ggplot2","ggpubr"))
library(readr)    
library(tidyr)    
library(dplyr)    
library(ggplot2)  
library(ggpubr)   

# 1) Read TPM matrix
tpm_df <- read_csv("NormalizedGeneExpression_TPM.csv",
                   col_types = cols(
                     ensembl_gene_id = col_character(),
                     Hugo_Symbol     = col_character(),
                     .default        = col_double()
                   )
)

# 2) Read Spearman-summary
corr_df <- read_csv("vegf_TPM.csv",
                    col_types = cols(
                      GeneID                           = col_character(),
                      GeneName                         = col_character(),
                      Spearman_Correlation_Coefficient = col_double(),
                      P_value                          = col_double(),
                      qvalue                           = col_double()
                    )
)

# 3) Read phenotype (mean_modz)
pheno_df <- read_csv("vegf_pheno.csv",
                     col_types = cols(
                       cytokine     = col_character(),
                       treatment    = col_character(),
                       donor_1_modz = col_double(),
                       donor_2_modz = col_double(),
                       donor_3_modz = col_double(),
                       mean_modz    = col_double()
                     )
)

# 4) Specify the two genes
genes_of_interest <- c(
  "ENSCAFG00000015332",  # MVB12A
  "ENSCAFG00000024217"   # TRAPPC5
)

# Manual gene name mapping
gene_name_map <- c(
  "ENSCAFG00000015332" = "MVB12A",
  "ENSCAFG00000024217" = "ENSCAFG00000024217 (TRAPPC5)"
)

# 5) Merge TPM + corr + keep only those two genes
plot_df <- tpm_df %>%
  inner_join(corr_df, by = c("ensembl_gene_id" = "GeneID")) %>%
  filter(ensembl_gene_id %in% genes_of_interest) %>%
  mutate(GeneName = gene_name_map[ensembl_gene_id]) %>%
  pivot_longer(
    cols = -c(
      ensembl_gene_id,
      Hugo_Symbol,
      GeneName,
      Spearman_Correlation_Coefficient,
      P_value,
      qvalue
    ),
    names_to  = "Sample",
    values_to = "TPM"
  ) %>%
  left_join(pheno_df, by = c("Sample" = "treatment")) %>%
  filter(!is.na(TPM) & !is.na(mean_modz)) %>%
  # Set factor order: MVB12A first, TRAPPC5 second
  mutate(GeneName = factor(GeneName, 
                           levels = c("MVB12A", "ENSCAFG00000024217 (TRAPPC5)")))

# 6) Build stats table for annotation
stats_df <- corr_df %>%
  filter(GeneID %in% genes_of_interest) %>%
  mutate(GeneName = gene_name_map[GeneID]) %>%
  dplyr::select(GeneName, rho = Spearman_Correlation_Coefficient,
                pval = P_value, qval = qvalue) %>%
  mutate(
    label = paste0(
      "ρ = ", round(rho,3), "\n",
      "q = ", signif(qval,3)
    ),
    # Set same factor order
    GeneName = factor(GeneName, 
                      levels = c("MVB12A", "ENSCAFG00000024217 (TRAPPC5)"))
  )

# 7) Plot with PLOS requirements
p <- ggplot(plot_df, aes(x = TPM, y = mean_modz)) +
  geom_point(size = 3, alpha = 0.8, color = "#2C3E50") +
  geom_smooth(method = "lm", se = TRUE,
              linetype = "dashed", color = "#E74C3C") +
  geom_text(
    data    = stats_df,
    aes(x = -Inf, y = -Inf, label = label),
    hjust   = 1.1, vjust = 1.1,
    size    = 3.5,
    color   = "black",
    family  = "Arial"
  ) +
  facet_wrap(~ GeneName, scales = "free") +
  labs(
    x = "Expression (TPM)",
    y = "Mean VEGF Modified Z-score"
  ) +
  theme_minimal(base_size = 12, base_family = "Arial") +
  theme(
    strip.text      = element_text(face = "bold.italic", size = 12),
    axis.title      = element_text(face = "bold"),
    axis.text       = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

# PLOS-compliant TIFF with LZW compression
ggsave(
  filename = "spearman_two_genes.tif",
  plot     = p,
  width    = 8,
  height   = 4,
  units    = "in",
  dpi      = 300,
  compression = "lzw",
  bg       = "white"
)

# EPS alternative
ggsave(
  filename = "spearman_two_genes.eps",
  plot     = p,
  width    = 8,
  height   = 4,
  units    = "in",
  dpi      = 300,
  device   = cairo_ps,
  bg       = "white"
)

print(p)

# ========================================
# EXTRACT LINEAR REGRESSION STATISTICS (for PLOS)
# ========================================

library(broom)

# Function to get complete regression stats
get_regression_stats <- function(gene_id, gene_name, data) {
  # Filter data for this gene
  gene_data <- data %>% 
    filter(ensembl_gene_id == gene_id)
  
  # Fit linear model
  lm_model <- lm(mean_modz ~ TPM, data = gene_data)
  
  # Extract coefficients with CIs
  coef_table <- broom::tidy(lm_model, conf.int = TRUE) %>%
    mutate(gene_id = gene_id, gene_name = gene_name) %>%
    select(gene_id, gene_name, term, estimate, std.error, statistic, p.value, conf.low, conf.high)
  
  # Extract model fit statistics
  fit_stats <- broom::glance(lm_model) %>%
    mutate(gene_id = gene_id, gene_name = gene_name) %>%
    select(gene_id, gene_name, r.squared, adj.r.squared, sigma, df, df.residual)
  
  return(list(coefficients = coef_table, fit = fit_stats))
}

# Get stats for both genes
mvb12a_stats <- get_regression_stats("ENSCAFG00000015332", "MVB12A", plot_df)
trappc5_stats <- get_regression_stats("ENSCAFG00000024217", "TRAPPC5", plot_df)

# Combine into publication tables
regression_coefficients <- bind_rows(
  mvb12a_stats$coefficients,
  trappc5_stats$coefficients
)

regression_fit <- bind_rows(
  mvb12a_stats$fit,
  trappc5_stats$fit
)

# Format for publication
regression_table <- regression_coefficients %>%
  left_join(regression_fit, by = c("gene_id", "gene_name")) %>%
  mutate(
    Gene = gene_name,
    Term = case_when(
      term == "(Intercept)" ~ "Intercept",
      term == "TPM" ~ "Slope (TPM)"
    ),
    Estimate = sprintf("%.4f", estimate),
    SE = sprintf("%.4f", std.error),
    `95% CI` = sprintf("[%.4f, %.4f]", conf.low, conf.high),
    t = sprintf("%.2f", statistic),
    `p-value` = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
    `R²` = sprintf("%.3f", r.squared),
    `Adj. R²` = sprintf("%.3f", adj.r.squared)
  ) %>%
  select(Gene, Term, Estimate, SE, `95% CI`, t, `p-value`, `R²`, `Adj. R²`) %>%
  # Remove duplicate R² rows
  mutate(
    `R²` = ifelse(Term == "Intercept", `R²`, ""),
    `Adj. R²` = ifelse(Term == "Intercept", `Adj. R²`, "")
  )

# Print and save
print(regression_table)
write_csv(regression_table, "S_Table_Fig3_Linear_Regression.csv")

cat("\n=== REGRESSION SUMMARY ===\n")
cat("MVB12A: R² =", sprintf("%.3f", mvb12a_stats$fit$r.squared), "\n")
cat("TRAPPC5: R² =", sprintf("%.3f", trappc5_stats$fit$r.squared), "\n")
