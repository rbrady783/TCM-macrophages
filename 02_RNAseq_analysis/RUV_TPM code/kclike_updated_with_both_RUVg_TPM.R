# Author:  Sunetra Das / Rachel Brady
# Date:   2024-02-28
#
# Project Name: kc_like_RUV_v2 with new dogs 2 and 3 - IMPROVED
#########################################################
#########################################################  
# Section: Reading in the input data
#########################################################
#########################################################  
# Getting the datasheet
gd <- read.csv("kc_like_modz.csv")
gd$Cell.Line[gd$Cell.Line %in% "D17"] <- "D-17"
gd$Cell.Line[gd$Cell.Line %in% "Den"] <- "DEN-HSA"
gd$Cell.Line[gd$Cell.Line %in% "OS2.4"] <- "OS2-4"
gd$Cell.Line[gd$Cell.Line %in%  "CML-6m"] <- "CML-6M"
# Moving the cell line column as rownames
rownames(gd) <- gd$Cell.Line
gd$Cell.Line <- NULL

# RUVg normalized data
expr <- read.csv("STAR_GeneCounts_RNAseq_UQ_RUVg.csv", check.names = F)
# Gene names
genenames <- subset(expr, select=c(1, 2))
genenames[1,]
# Changes the cell line name syntax to match 
rownames(expr) <- expr$ensembl_gene_id
expr$ensembl_gene_id <- NULL
expr$Hugo_Symbol <- NULL
expr[c(1:10),]

#########################################################
# Section: Data Quality Checks (NEW)
#########################################################
# Verify sample matching
cat("Cytokine samples:", paste(rownames(gd), collapse = ", "), "\n")
cat("RNA-seq samples:", paste(colnames(expr), collapse = ", "), "\n")

# Check for mismatches
missing_in_expr <- setdiff(rownames(gd), colnames(expr))
missing_in_cytokine <- setdiff(colnames(expr), rownames(gd))
if(length(missing_in_expr) > 0 | length(missing_in_cytokine) > 0) {
  warning("Sample mismatch detected!")
  print(missing_in_expr)
  print(missing_in_cytokine)
}

# Filter out zero variance genes
gene_vars <- apply(expr, 1, function(x) var(x, na.rm = TRUE))
zero_var_genes <- sum(gene_vars == 0 | is.na(gene_vars))
cat("Genes with zero variance:", zero_var_genes, "\n")
valid_genes <- gene_vars > 0 & !is.na(gene_vars)
expr <- expr[valid_genes, ]
genenames <- genenames[genenames$ensembl_gene_id %in% rownames(expr), ]

#########################################################
#########################################################
# Section: Running the correlations in loops
#########################################################
######################################################### 
# Creating empty data set for logging results
df <- data.frame("GeneID" = genenames$ensembl_gene_id,
                 "GeneName"=genenames$Hugo_Symbol, 
                 "Spearman_Correlation_Coefficient" = NA, 
                 "P_value"=NA, stringsAsFactors = F)

# Analyte data
a <- subset(gd, select=4)

# FIXED: Use correct number of genes
for (i in 1:nrow(expr)){
  
  # Correlation function
  # Gene name and symbol
  ge <- rownames(expr)[i]
  gn <- genenames$Hugo_Symbol[genenames$ensembl_gene_id %in% ge]
  
  # IMPROVED: Add error handling
  tryCatch({
    # Gene expression
    b <- data.frame(t(expr[rownames(expr) %in% ge,]))
    # Merging the analyte data with gene expression 
    c <- merge(a, b, by=0)
    
    # Check for sufficient data
    if(nrow(c) < 3) {
      warning(paste("Insufficient data for gene", ge))
      next
    }
    
    # Correlation test
    d <- cor.test(c[,2], c[,3], method = "spearman")
    # FIXED: Use result from cor.test instead of calculating twice
    df$Spearman_Correlation_Coefficient[df$GeneID %in% ge] <- d$estimate
    df$P_value[df$GeneID %in% ge] <- d$p.value
    
  }, error = function(e) {
    warning(paste("Correlation failed for gene", ge, ":", e$message))
  })
  
  # print data
  print(i)
  print(df[i,])
  
}  

# Deleting the genes where the correlation could not be calculated
df <- df[!is.na(df$P_value),]

#########################################################
# Section: FDR Correction - BOTH METHODS (NEW)
#########################################################

# Method 1: Your original fdrtool method
fdis <- fdrtool::fdrtool(df$P_value, statistic=c("pvalue"),
                         plot=TRUE, color.figure=TRUE, 
                         verbose=TRUE, cutoff.method=c("pct0"), pct0=0.75)

# Method 2: Benjamini-Hochberg correction
bh_qvalues <- p.adjust(df$P_value, method = "BH")

# Binding both qvalue methods to main data
fcor.sig <- cbind(df, fdis$qval, bh_qvalues)
# renaming headers
colnames(fcor.sig)[length(colnames(fcor.sig))-1] <- "fdrtool_qvalue"
colnames(fcor.sig)[length(colnames(fcor.sig))] <- "BH_qvalue"

# Ordering the data
fcor.sig <- fcor.sig[order(fcor.sig$P_value),]
head(fcor.sig)

# Results summary
length(fcor.sig$GeneID[fcor.sig$P_value < 0.05])
length(fcor.sig$GeneID[fcor.sig$fdrtool_qvalue < 0.05])
length(fcor.sig$GeneID[fcor.sig$BH_qvalue < 0.05])
range(fcor.sig$P_value)

# Saving the data
fn <- paste0("C:/Users/brady/OneDrive/Desktop/TCM_RNA/Output/CorrelationTables/kc_like_RUVg_improved.csv")
write.csv(fcor.sig, fn, row.names = FALSE)

#########################################################
#########################################################
#########################################################

# Date:   2024-04-07
#
# Project Name: kc_like_TPM - IMPROVED
#########################################################
#########################################################  
# Section: Reading in the input data
#########################################################
#########################################################  
# Getting the datasheet
gd <- read.csv("kc_like_modz.csv")
gd$Cell.Line[gd$Cell.Line %in% "D17"] <- "D-17"
gd$Cell.Line[gd$Cell.Line %in% "Den"] <- "DEN-HSA"
gd$Cell.Line[gd$Cell.Line %in% "OS2.4"] <- "OS2-4"
gd$Cell.Line[gd$Cell.Line %in%  "CML-6m"] <- "CML-6M"
# Moving the cell line column as rownames
rownames(gd) <- gd$Cell.Line
gd$Cell.Line <- NULL

# TPM normalized data
expr <- read.csv("NormalizedGeneExpression_TPM.csv", check.names = F)
# Gene names
genenames <- subset(expr, select=c(1, 2))
genenames[1,]
# Changes the cell line name syntax to match 
rownames(expr) <- expr$ensembl_gene_id
expr$ensembl_gene_id <- NULL
expr$Hugo_Symbol <- NULL
expr[c(1:10),]

#########################################################
# Section: Data Quality Checks (NEW)
#########################################################
# Verify sample matching
cat("Cytokine samples:", paste(rownames(gd), collapse = ", "), "\n")
cat("RNA-seq samples:", paste(colnames(expr), collapse = ", "), "\n")

# Check for mismatches
missing_in_expr <- setdiff(rownames(gd), colnames(expr))
missing_in_cytokine <- setdiff(colnames(expr), rownames(gd))
if(length(missing_in_expr) > 0 | length(missing_in_cytokine) > 0) {
  warning("Sample mismatch detected!")
  print(missing_in_expr)
  print(missing_in_cytokine)
}

# Filter out zero variance genes
gene_vars <- apply(expr, 1, function(x) var(x, na.rm = TRUE))
zero_var_genes <- sum(gene_vars == 0 | is.na(gene_vars))
cat("Genes with zero variance:", zero_var_genes, "\n")
valid_genes <- gene_vars > 0 & !is.na(gene_vars)
expr <- expr[valid_genes, ]
genenames <- genenames[genenames$ensembl_gene_id %in% rownames(expr), ]

#########################################################
#########################################################
# Section: Running the correlations in loops
#########################################################
######################################################### 
# Creating empty data set for logging results
df <- data.frame("GeneID" = genenames$ensembl_gene_id,
                 "GeneName"=genenames$Hugo_Symbol, 
                 "Spearman_Correlation_Coefficient" = NA, 
                 "P_value"=NA, stringsAsFactors = F)

# Analyte data
a <- subset(gd, select=4)

# FIXED: Use correct number of genes
for (i in 1:nrow(expr)){
  
  # Correlation function
  # Gene name and symbol
  ge <- rownames(expr)[i]
  gn <- genenames$Hugo_Symbol[genenames$ensembl_gene_id %in% ge]
  
  # IMPROVED: Add error handling
  tryCatch({
    # Gene expression
    b <- data.frame(t(expr[rownames(expr) %in% ge,]))
    # Merging the analyte data with gene expression 
    c <- merge(a, b, by=0)
    
    # Check for sufficient data
    if(nrow(c) < 3) {
      warning(paste("Insufficient data for gene", ge))
      next
    }
    
    # Correlation test
    d <- cor.test(c[,2], c[,3], method = "spearman")
    # FIXED: Use result from cor.test instead of calculating twice
    df$Spearman_Correlation_Coefficient[df$GeneID %in% ge] <- d$estimate
    df$P_value[df$GeneID %in% ge] <- d$p.value
    
  }, error = function(e) {
    warning(paste("Correlation failed for gene", ge, ":", e$message))
  })
  
  # print data
  print(i)
  print(df[i,])
  
}  

# Deleting the genes where the correlation could not be calculated
df <- df[!is.na(df$P_value),]

#########################################################
# Section: FDR Correction - BOTH METHODS (NEW)
#########################################################

# Method 1: Your original fdrtool method
fdis <- fdrtool::fdrtool(df$P_value, statistic=c("pvalue"),
                         plot=TRUE, color.figure=TRUE, 
                         verbose=TRUE, cutoff.method=c("pct0"), pct0=0.75)

# Method 2: Benjamini-Hochberg correction
bh_qvalues <- p.adjust(df$P_value, method = "BH")

# Binding both qvalue methods to main data
fcor.sig <- cbind(df, fdis$qval, bh_qvalues)
# renaming headers
colnames(fcor.sig)[length(colnames(fcor.sig))-1] <- "fdrtool_qvalue"
colnames(fcor.sig)[length(colnames(fcor.sig))] <- "BH_qvalue"

# Ordering the data
fcor.sig <- fcor.sig[order(fcor.sig$P_value),]
head(fcor.sig)

# Results summary
length(fcor.sig$GeneID[fcor.sig$P_value < 0.05])
length(fcor.sig$GeneID[fcor.sig$fdrtool_qvalue < 0.05])
length(fcor.sig$GeneID[fcor.sig$BH_qvalue < 0.05])
range(fcor.sig$P_value)

# Saving the data
fn <- paste0("C:/Users/brady/OneDrive/Desktop/TCM_RNA/Output/CorrelationTables/kc_like_TPM_improved.csv")
write.csv(fcor.sig, fn, row.names = FALSE)
