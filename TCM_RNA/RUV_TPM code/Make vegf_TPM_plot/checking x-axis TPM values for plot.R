library(readr)
library(dplyr)
library(tidyr)

# Read the TPM data
tpm_df <- read_csv("NormalizedGeneExpression_TPM.csv",
                   col_types = cols(
                     ensembl_gene_id = col_character(),
                     Hugo_Symbol     = col_character(),
                     .default        = col_double()
                   )
)

# Extract the two genes
genes_of_interest <- c("ENSCAFG00000015332", "ENSCAFG00000024217")

# Filter and reshape to see the TPM values
tpm_extract <- tpm_df %>%
  filter(ensembl_gene_id %in% genes_of_interest) %>%
  pivot_longer(
    cols = -c(ensembl_gene_id, Hugo_Symbol),
    names_to = "Sample",
    values_to = "TPM"
  ) %>%
  arrange(ensembl_gene_id, TPM)

# Print results for MVB12A (ENSCAFG00000015332)
cat("\n=== MVB12A (ENSCAFG00000015332) TPM values (lowest to highest) ===\n")
mvb12a <- tpm_extract %>% filter(ensembl_gene_id == "ENSCAFG00000015332")
print(mvb12a %>% select(Sample, TPM), n = 100)
cat("\nRange:", min(mvb12a$TPM), "to", max(mvb12a$TPM), "\n")

# Print results for second gene (ENSCAFG00000024217)
cat("\n=== ENSCAFG00000024217 (TRAPPC5) TPM values (lowest to highest) ===\n")
gene2 <- tpm_extract %>% filter(ensembl_gene_id == "ENSCAFG00000024217")
print(gene2 %>% select(Sample, TPM), n = 100)
cat("\nRange:", min(gene2$TPM), "to", max(gene2$TPM), "\n")

# Summary statistics
cat("\n=== Summary Statistics ===\n")
summary_stats <- tpm_extract %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    Min = min(TPM),
    Q1 = quantile(TPM, 0.25),
    Median = median(TPM),
    Mean = mean(TPM),
    Q3 = quantile(TPM, 0.75),
    Max = max(TPM),
    N = n()
  )
print(summary_stats)