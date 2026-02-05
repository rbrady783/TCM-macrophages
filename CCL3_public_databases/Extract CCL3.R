# Load required libraries
library(dplyr)
library(readr)
library(ggplot2)

# Get a list of all directories in the current working directory
all_dirs <- list.dirs(recursive = FALSE, full.names = TRUE)

# Filter to include only those ending with "_quant"
quant_folders <- all_dirs[grep("_quant$", all_dirs)]

# Check what folders we have:
print(quant_folders)

# Function to extract CCL3 info from a quant.sf file in a given folder
extract_CCL3 <- function(folder) {
  quant_file <- file.path(folder, "quant.sf")
  if (!file.exists(quant_file)) {
    warning(paste("File not found:", quant_file))
    return(NULL)
  }
  # Read the quant file (tab-delimited)
  df <- read_tsv(quant_file)
  
  # Debugging: Print unique names to check transcript IDs
  cat("Unique names in", folder, ":\n")
  print(unique(df$Name))
  
  # Filter for rows corresponding to CCL3 using grepl for flexibility
  df_CCL3 <- df %>% filter(grepl("ENSCAFT00000028854", Name))
  
  # Add a column for sample name (extracted from folder name)
  sample_name <- basename(folder)
  df_CCL3 <- df_CCL3 %>% mutate(Sample = sample_name)
  
  return(df_CCL3)
}

# Loop over the filtered quant folders and bind results together
ccl3_list <- lapply(quant_folders, extract_CCL3)
ccl3_data <- bind_rows(ccl3_list)

# Save the extracted data to a CSV file
write_csv(ccl3_data, "CCL3_expression.csv")

# Create a bar plot of TPM levels for each sample
ggplot(ccl3_data, aes(x = Sample, y = TPM)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(title = "CCL3 Expression",
       x = "Sample",
       y = "TPM") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
