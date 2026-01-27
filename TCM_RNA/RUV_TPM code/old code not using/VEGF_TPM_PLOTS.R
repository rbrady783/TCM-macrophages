Project Name VEGF core
#
# Script Description: Correlation plots protein VEGF conc with canine transcriptome

#########################################################
#########################################################  


#########################################################
#########################################################  

#########################################################
#########################################################

# Section: Libraries and functions to load

#########################################################
#########################################################  
library(stringi)
library(ggpp)
library(ggpmisc)
library(ggrepel)
library(gridExtra)
library(grid)
library(fdrtool)
library(colorBlindness)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggdendro)
library(ggpubr)
library(scales)
options(stringsAsFactors = FALSE)
options(scipen=999)
date <- Sys.Date()
date <- gsub("-", ".", date)


#My theme for ggplots
my_Theme = theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 12, angle =0,
                             hjust = 1, vjust = 0.5),
  axis.title.y = element_text(size = 12),
  axis.text.y = element_text(size=12),# 12
  legend.title=element_text(size=10),
  legend.text=element_text(size=10),
  strip.text = element_text(size=13),
  panel.grid.major.y = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.background = element_rect(fill = NA),
  axis.line = element_line(color = "black"),
  axis.text = element_text(colour = "black"),
  strip.background=element_rect(color = "black",  fill="grey75")
)
#########################################################
#########################################################

# Section: Reading in data

#########################################################
#########################################################  

# Getting the VEGF
gd <- xlsx::read.xlsx("VEGF to send to Sunetra.xlsx", sheetIndex = 1)
gd$Cell.Line[gd$Cell.Line %in% "D17"] <- "D-17"
gd$Cell.Line[gd$Cell.Line %in% "Den"] <- "DEN-HSA"
gd$Cell.Line[gd$Cell.Line %in% "OS2.4"] <- "OS2-4"
gd$Cell.Line[gd$Cell.Line %in%  "CML-6m"] <- "CML-6M"
# Moving the cell line column as rownames
rownames(gd) <- gd$Cell.Line
gd$Cell.Line <- NULL

# RUGg normalized data
expr <- read.csv("NormalizedGeneExpression_23_Canine_CellLines_TPM.csv", check.names = F)

# Gene names
genenames <- subset(expr, select=c(1, 2))
genenames[1,]

# Changes the cell line name syntax to match with Dan's format
rownames(expr) <- expr$ensembl_gene_id
expr$ensembl_gene_id <- NULL
expr$Hugo_Symbol <- NULL
expr[c(1:10),]

#########################################################
#########################################################

# Section: Plotting all data in a loop

#########################################################
#########################################################  
# VEGF data
a <- subset(gd, select=4)

#Plot empty list
p <- list()

# Loop
for (j in 1:length(expr$CLL1390)){
  # Gene name and symbol
  ge <- rownames(expr)[j]
  gn <- genenames$Hugo_Symbol[genenames$ensembl_gene_id %in% ge]
  # Gene expression
  b <- data.frame(t(expr[rownames(expr) %in% ge,]))
  
  # Plotting the VEGF data vs gene expression
  a$TSG <- rownames(a)
  firstP <- merge(a, b, by.y=0, by.x="TSG" )
  colnames(firstP)[2:3] <- c("VEGF", "GeneName")
  
  #Filename
  fn <- paste0("C:/Users/brady/OneDrive/Desktop/TCM_RNA/Output/CorrelationPlots", gn, "_", ge, ".pdf")
  # Saving plot
  p[[j]] <-  ggplot(firstP, aes(x=VEGF, y=GeneName)) +
    geom_label_repel(aes(label = TSG),
                     box.padding   = 0.35,
                     point.padding = 0.5,
                     segment.color = 'red', 
                     size = 1) +
    geom_point() +
    xlab("VEGF")+
    ylab(paste0("RUVg normalized \n gene expression")) +
    ggtitle(paste0(gn, "_", ge))+
    theme(plot.title = element_text(color="blue", size=10, face="bold"))+
    stat_poly_line() +
    stat_poly_eq(vjust = 3,use_label(c("R2", "P")), size=2) +
    scale_color_manual(values=c( "purple4" , "black"))+
    my_Theme +
    theme(legend.position="none")
  ggsave(p[[j]], file=fn, width = 6, height = 4, units = "in")
  print(j)
  print(gn)
}
