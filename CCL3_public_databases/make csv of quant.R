
library(readr)
library(ggplot2)

# Replace the path with the appropriate folder and file name
df <- read_tsv("DRR532498_quant/quant.sf")

# Write the data frame to a CSV file
write_csv(df, "DRR532498_quant/quant.csv")

#plot
# Read the CSV file containing your expression data
ccl3_data <- read_csv("CCL3_expression.csv")

# Print the data to inspect it
print(ccl3_data)

library(dplyr)

# Recode sample names (adjust the mappings to your desired names)
ccl3_data <- ccl3_data %>%
  mutate(Sample = recode(Sample,
                         "DRR532498_quant" = "CHS1",
                         "DRR532499_quant" = "CHS2",
                         "DRR532500_quant" = "CHS3",
                         "DRR532501_quant" = "CHS4",
                         "DRR532502_quant" = "CHS5",
                         "DRR532503_quant" = "CHS6",
                         "DRR532504_quant" = "CHS7",
                         "DRR532505_quant" = "CHS8",
                         "DRR532506_quant" = "MHT2",
                         "DRR532507_quant" = "DHS1",
                         "DRR532508_quant" = "DHS2",
                         "DRR532509_quant" = "DH82"
  ))

# Create a bar plot of TPM levels for each sample
ggplot(ccl3_data, aes(x = Sample, y = TPM)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "CCL3 Expression Levels (TPM) Across Samples",
       x = "Sample",
       y = "TPM")

