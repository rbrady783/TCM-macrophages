#— Install once if needed:
# install.packages(c("readr","tidyr","dplyr","ggplot2","ggpubr"))

library(readr)    # read_csv()
library(tidyr)    # pivot_longer()
library(dplyr)    # data wrangling
library(ggplot2)  # plotting
library(ggpubr)   # stat_cor()

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
  "ENSCAFG00000024217"   # second gene
)

# 5) Merge TPM + corr + keep only those two genes
plot_df <- tpm_df %>%
  inner_join(corr_df, by = c("ensembl_gene_id" = "GeneID")) %>%
  filter(ensembl_gene_id %in% genes_of_interest) %>%
  # reshape TPM → long form:
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
  # join in phenotype
  left_join(pheno_df, by = c("Sample" = "treatment")) %>%
  filter(!is.na(TPM) & !is.na(mean_modz))

# 6) Build a small stats table for annotation
stats_df <- corr_df %>%
  filter(GeneID %in% genes_of_interest) %>%
  dplyr::select(GeneName, rho = Spearman_Correlation_Coefficient,
         pval = P_value, qval = qvalue) %>%
  mutate(
    label = paste0(
      "ρ = ", round(rho,3), "\n",
      "q = ", signif(qval,3)
    )
  )

# 7) Plot with facet_wrap 
p <- ggplot(plot_df, aes(x = TPM, y = mean_modz)) +
  geom_point(size = 3, alpha = 0.8, color = "#2C3E50") +
  geom_smooth(method = "lm", se = TRUE,
              linetype = "dashed", color = "#E74C3C") +
  # add per-panel text from stats_df
  geom_text(
    data    = stats_df,
    aes(x = -Inf, y = -Inf, label = label),
    hjust   = 1.1, vjust = 1.1,
    size    = 4,
    color   = "black"
  ) +
  facet_wrap(~ GeneName, scales = "free") +
  labs(
    x = "Expression (TPM)",
    y = "Mean VEGF Modified Z-score"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text      = element_text(face = "bold", size = 14),
    axis.title      = element_text(face = "bold"),
    axis.text       = element_text(color = "black"),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = "spearman_two_genes.png",
  plot     = p,
  width    = 8,       # in inches
  height   = 4,       # in inches
  units    = "in",
  dpi      = 600,
  bg       = "transparent"
)

ggsave(
  filename = "spearman_two_genes.eps",
  plot     = p,
  width    = 8,       # in inches
  height   = 4,       # in inches
  units    = "in",
  dpi      = 300,
  device   = cairo_ps,
  bg       = "white"
)

print(p)
