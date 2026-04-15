# Load required libraries
library(readr)
library(dplyr)
library(openxlsx)

# Set your working directory to the GSEA folder
# setwd("path/to/your/GSEA/folder")

# Define the cytokine mapping based on your file names
cytokine_files <- list(
  "CCL2" = "CCL2_GSEA_All_Collections.csv",
  "IL8" = "il.8_GSEA_All_Collections.csv", 
  "IL10" = "il.10_GSEA_All_Collections.csv",
  "KC_like" = "kc.like_GSEA_All_Collections.csv",
  "TGFB" = "tgf.b_GSEA_All_Collections.csv",
  "TNFA" = "tnf.a_GSEA_All_Collections.csv",
  "VEGF" = "vegf_GSEA_All_Collections.csv"
)

# Function to read and clean GSEA data
read_gsea_data <- function(file_path, cytokine_name) {
  tryCatch({
    data <- read_csv(file_path, show_col_types = FALSE)
    
    # Add cytokine column for identification
    data$Cytokine <- cytokine_name
    
    # Clean column names (remove spaces, special characters)
    colnames(data) <- gsub("[^A-Za-z0-9_]", "_", colnames(data))
    colnames(data) <- gsub("_+", "_", colnames(data))
    colnames(data) <- gsub("^_|_$", "", colnames(data))
    
    return(data)
  }, error = function(e) {
    warning(paste("Could not read file:", file_path, "- Error:", e$message))
    return(NULL)
  })
}

# Read all GSEA files
gsea_data_list <- list()
combined_data <- data.frame()

cat("Reading GSEA files...\n")
for(cytokine in names(cytokine_files)) {
  file_path <- cytokine_files[[cytokine]]
  
  if(file.exists(file_path)) {
    cat(paste("Reading:", cytokine, "from", file_path, "\n"))
    
    # Read the data
    data <- read_gsea_data(file_path, cytokine)
    
    if(!is.null(data)) {
      gsea_data_list[[cytokine]] <- data
      
      # Add to combined dataset
      if(nrow(combined_data) == 0) {
        combined_data <- data
      } else {
        # Combine datasets, handling different column structures
        combined_data <- bind_rows(combined_data, data)
      }
    }
  } else {
    warning(paste("File not found:", file_path))
  }
}

# Display summary
cat("\nSummary of data read:\n")
for(cytokine in names(gsea_data_list)) {
  cat(paste(cytokine, ":", nrow(gsea_data_list[[cytokine]]), "pathways\n"))
}

# Create Excel workbook with separate sheets for each cytokine
wb <- createWorkbook()

# Add each cytokine as a separate sheet
for(cytokine in names(gsea_data_list)) {
  addWorksheet(wb, cytokine)
  writeData(wb, cytokine, gsea_data_list[[cytokine]])
  
  # Auto-adjust column widths
  setColWidths(wb, cytokine, cols = 1:ncol(gsea_data_list[[cytokine]]), widths = "auto")
}

# Save the Excel file
excel_filename <- paste0("Combined_GSEA_Results_", Sys.Date(), ".xlsx")
saveWorkbook(wb, excel_filename, overwrite = TRUE)
