#########################################################
#########################################################
#########################################################

# Author:  Sunetra Das / Rachel Brady
# # Date:   2024-02-28
#
# Project Name: VEGF cor
#
# Script Description: Correlation protein VEGF conc with canine transcriptome

#########################################################
#########################################################  


#########################################################
#########################################################

# Section: Reading in the input data

#########################################################
#########################################################  


# Getting the VEGF datasheet
gd <- read.csv("vegf_modz.csv")
gd$Cell.Line[gd$Cell.Line %in% "D17"] <- "D-17"
gd$Cell.Line[gd$Cell.Line %in% "Den"] <- "DEN-HSA"
gd$Cell.Line[gd$Cell.Line %in% "OS2.4"] <- "OS2-4"
gd$Cell.Line[gd$Cell.Line %in%  "CML-6m"] <- "CML-6M"

# Moving the cell line column as rownames
rownames(gd) <- gd$Cell.Line
gd$Cell.Line <- NULL

# RUGg normalized data
expr <- read.csv("STAR_GeneCounts_RNAseq_UQ_RUVg.csv", check.names = F)

# Gene names
genenames <- subset(expr, select=c(1, 2))
genenames[1,]

# Changes the cell line name syntax to match 
rownames(expr) <- expr$ensembl_gene_id
expr$ensembl_gene_id <- NULL
expr$Hugo_Symbol <- NULL
expr[c(1:10),]


# Section: Spearman correlations : Growth rate vs gene expression

#########################################################
######################################################### 

# Creating empty data set for logging results
df <- data.frame("GeneID" = genenames$ensembl_gene_id,
                 "GeneName"=genenames$Hugo_Symbol, 
                 "Spearman_Correlation_Coefficient" = NA, 
                 "P_value"=NA, stringsAsFactors = F)

# VEGF data
a <- subset(gd, select=4)

for (i in 1:length(expr$Abrams)){
  
  # Correlation function
  # Gene name and symbol
  ge <- rownames(expr)[i]
  gn <- genenames$Hugo_Symbol[genenames$ensembl_gene_id %in% ge]
  # Gene expression
  b <- data.frame(t(expr[rownames(expr) %in% ge,]))
  # Merging the VEGF data with gene expression 
  c <- merge(a, b, by=0)
  # Correlation test
  d <- cor.test(c[,2], c[,3], method = "spearman")
  # Adding the correlation data to output data set
  df$Spearman_Correlation_Coefficient[df$GeneID %in% ge] <- cor(c[,2], c[,3], method = "spearman")
  df$P_value[df$GeneID %in% ge] <- d$p.value
  # print data
  print(i)
  print(df[i,])
  
}  

# Deleting the genes where the correlation could not be calculated
df <- df[!is.na(df$P_value),]

# Adding FDR
fdis <- fdrtool::fdrtool(df$P_value, statistic=c("pvalue"),
                         plot=TRUE, color.figure=TRUE, 
                         verbose=TRUE, cutoff.method=c("pct0"), pct0=0.75)

# Binding the qvalue to main data
fcor.sig <- cbind(df,  fdis$qval)

# renaming header
colnames(fcor.sig)[length(colnames(fcor.sig))] <- "qvalue"

# Ordering the data
fcor.sig <- fcor.sig[order(fcor.sig$P_value),]
head(fcor.sig)

length(fcor.sig$GeneID[fcor.sig$P_value <0.05])
range(fcor.sig$P_value)

# Saving the data
fn <- paste0("C:/Users/brady/OneDrive/Desktop/TCM_RNA/Output/CorrelationTables/Spearman_Correlation_RUVg_with_VEGF_modz.csv")
write.csv(fcor.sig, fn, row.names = FALSE)


