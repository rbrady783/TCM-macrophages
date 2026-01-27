# Load required libraries
library(readr)
library(dplyr)
library(openxlsx)
library(stringr)

#install.packages("openxlsx")
# Set working directory (adjust path as needed)
# setwd("C:/path/to/your/CorrelationTables/folder")

# Define function to clean file names
clean_file_name <- function(filename) {
  # Remove "_improved" from the name
  cleaned <- str_replace(filename, "_improved", "")
  # Remove file extensions
  cleaned <- str_replace(cleaned, "\\.(csv|xlsx?)$", "")
  return(cleaned)
}

# Get list of all CSV files with _RUVg or _TPM in the name
csv_files <- list.files(pattern = "(_RUVg|_TPM).*\\.(csv|xlsx?)$", full.names = TRUE)

# Print found files
cat("Found", length(csv_files), "files:\n")
print(csv_files)

# Initialize list to store data frames
data_list <- list()

# Read each file and store with cleaned names
for (file in csv_files) {
  cat("Processing:", file, "\n")
  
  # Get just the filename without path
  file_name <- basename(file)
  
  # Clean the name
  clean_name <- clean_file_name(file_name)
  
  # Read the file (handles both CSV and Excel)
  if (str_detect(file, "\\.csv$")) {
    data <- read_csv(file, show_col_types = FALSE)
  } else {
    data <- read.xlsx(file, sheet = 1)
  }
  
  # Add source column to track which file data came from
  data$source_file <- clean_name
  
  # Store in list with cleaned name
  data_list[[clean_name]] <- data
  
  cat("  -> Cleaned name:", clean_name, "\n")
  cat("  -> Dimensions:", nrow(data), "rows x", ncol(data), "columns\n\n")
}

# Create a workbook with separate sheets for each file
wb <- createWorkbook()

# Add each dataset as a separate sheet
for (name in names(data_list)) {
  addWorksheet(wb, sheetName = name)
  writeData(wb, sheet = name, x = data_list[[name]])
}

# Also create a combined sheet if all files have similar structure
# First, let's check if we can combine them
cat("Checking if files can be combined...\n")

# Get column names from each file
col_names_list <- lapply(data_list, names)

# Check if all files have the same columns (excluding source_file)
all_cols <- unique(unlist(col_names_list))
common_cols <- Reduce(intersect, lapply(col_names_list, function(x) x[x != "source_file"]))

cat("Common columns across all files:", length(common_cols), "\n")
cat("Total unique columns:", length(all_cols), "\n")

# If files have similar structure, create combined sheet
if (length(common_cols) > 0) {
  # Combine all data, keeping only common columns
  combined_data <- bind_rows(lapply(data_list, function(df) {
    df %>% select(all_of(c(common_cols, "source_file")))
  }))
  
  # Add combined sheet
  addWorksheet(wb, sheetName = "Combined_All_Data")
  writeData(wb, sheet = "Combined_All_Data", x = combined_data)
  
  cat("Created combined sheet with", nrow(combined_data), "total rows\n")
}

# Save the Excel file
output_file <- "Combined_Correlation_Data.xlsx"
saveWorkbook(wb, output_file, overwrite = TRUE)

cat("\n=== SUMMARY ===\n")
cat("Successfully processed", length(data_list), "files\n")
cat("Output saved as:", output_file, "\n")
cat("Sheets created:\n")
for (name in names(data_list)) {
  cat("  -", name, "(", nrow(data_list[[name]]), "rows )\n")
}
if (exists("combined_data")) {
  cat("  - Combined_All_Data (", nrow(combined_data), "rows )\n")
}

# Display summary of name changes
cat("\n=== NAME CLEANING SUMMARY ===\n")
original_names <- basename(csv_files)
cleaned_names <- sapply(original_names, clean_file_name)

for (i in seq_along(original_names)) {
  cat(sprintf("%-25s -> %s\n", original_names[i], cleaned_names[i]))
}

cat("\nDone! Check the file:", output_file, "\n")
