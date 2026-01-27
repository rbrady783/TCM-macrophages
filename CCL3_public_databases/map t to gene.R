if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("tximport")

library(tximport)
library(readr)

# Assume you have a CSV mapping file with columns: transcript_id,gene_id
tx2gene <- read_csv("map_ccl3.csv") 

# List the Salmon quant.sf files (adjust path/pattern as needed)
salmon_files <- list.files("C:/Users/brady/OneDrive/Desktop/CCL3", 
                           pattern = "quant.sf", 
                           full.names = TRUE, 
                           recursive = TRUE)
names(salmon_files) <- basename(dirname(salmon_files))

# Check the list of files
print(salmon_files)

# Import Salmon quantifications and summarize to gene level
txi <- tximport(salmon_files, type = "salmon", tx2gene = tx2gene)
salmon_gene_counts <- txi$counts
salmon_tpm <- txi$abundance
