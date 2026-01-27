# Install and load AnnotationHub
# Install AnnotationHub from Bioconductor
#BiocManager::install("AnnotationHub")
#library(AnnotationHub)

ccl2_path <- dplyr::select(ccl2_pairwise, 2, 4, 8)

# Assuming your data frame is 'ccl2_pairwise' and ENSCAF IDs are in the 'names' column
enscaf_ids <- ccl2_path$Hugo_Symbol  # Extract ENSCAF gene IDs

# List available keytypes for org.Cf.eg.db
keytypes(org.Cf.eg.db)

mapped_ids <- select(org.Cf.eg.db, 
                     keys = enscaf_ids,       # ENSCAF gene IDs or symbols
                     columns = c("ENTREZID"), # Columns to fetch
                     keytype = "SYMBOL")      # Keytype for gene symbols

# Check the results
head(mapped_ids)
