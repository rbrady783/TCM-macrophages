library(tidyverse)

# Read the CSV file
data_wide <- read_csv("CCL3_facc.csv")

# Filter for your gene of interest and pivot the sample columns into long format.
gene_data <- data_wide %>%
  filter(ensembl_gene_id == "ENSCAFG00000018167") %>%
  pivot_longer(
    cols = -c(ensembl_gene_id, Hugo_Symbol),  # Only pivot the sample columns
    names_to = "Sample",
    values_to = "TPM"
  )

# View the result
print(gene_data)

ggplot(gene_data, aes(x = Sample, y = TPM)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Expression of CCL3 Across Samples",
       x = "Sample",
       y = "TPM")
