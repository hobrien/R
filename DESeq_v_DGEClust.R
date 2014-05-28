library(ggplot2)
library(plyr)
library(DESeq)
library(VennDiagram)

setwd("/Users/HeathOBrien/")

#read in raw count data
selaginellaCountTable =  read.table( file="Bioinformatics/Selaginella/Counts/Selag_counts.txt" , header=TRUE, row.names=1 )

#Create a set of factors to discribe the tissues
condition = factor( c( "leaf1", "leaf2", "leaf3", "leaf4" ) )

#subset the KRAUS data
KRAUS_CountTable = selaginellaCountTable [c('KRAUS1', 'KRAUS2', 'KRAUS3','KRAUS4')]

#create DESeq object from KRAUS count data and conditions and estimate dipersions
KRAUS_cds<-newCountDataSet( KRAUS_CountTable, condition)
KRAUS_cds <- estimateSizeFactors(KRAUS_cds)
KRAUS_cds <- estimateDispersions(KRAUS_cds, method="blind", sharingMode="fit-only")

# =======================================================================================
#test for differential expression of KRAUS1 and KRAUS2
DEseq<- nbinomTest(KRAUS_cds, "leaf1", "leaf2")

#rename columns and add info about species and comparison
DEseq["samples"]<-"12"
DEseq["species"]<-"KRAUS"
DEseq<-rename(DEseq, c("pval" = "DESeq_pval", "padj"="DESeq_padj"))

#read in results from DGEClust
DGEClust<-read.table(file="Bioinformatics/Selaginella/DGEClust/KRAUS12_pvals.txt", header=T)

#sort rows
DGEClust<-DGEClust[order(as.numeric(rownames(DGEClust))),]

#add DGESClust data to results
DEseq["DGEClust_pval"]<-DGEClust["pval"]
DEseq["DGEClust_padj"]<-DGEClust["padj"]

write.table(DEseq[DEseq$DESeq_padj < 0.1 & DEseq$DGEClust_padj < 0.1, ], 
  "Google Drive/Selaginella/DGEClust/KRAUS12_overlap.txt")

#this will replace extreme values with zero, which makes the plot nicer
DEseq[DEseq$DESeq_padj < 1e-5, ]["DESeq_padj"]<-0

#make scatterplot of adjusted p-values and Venn Diagram of overlap
pdf(file="Google Drive/Selaginella/DGEClust/KRAUS12.pdf")
ggplot(DEseq, aes(x=DESeq_padj, y=DGEClust_padj)) +
  geom_point(size=1.2) +
  scale_x_log10() +
  scale_y_log10() +
  geom_hline(yintercept=0.1, size=.2) +
  geom_vline(xintercept=0.1, size=.2) +
  theme_bw() +
  ggtitle("KRAUS12")
  
draw.pairwise.venn(area1=nrow(DEseq[DEseq$DGEClust_padj < 0.1, ]), 
  area2=nrow(DEseq[DEseq$DESeq_padj < 0.1, ]), 
  cross.area=nrow(DEseq[DEseq$DESeq_padj < 0.1 & DEseq$DGEClust_padj < 0.1, ]), 
  fill=c('blue', 'red'), 
  cat.col = c('blue', 'red'),
  category=c("DESeq", "DGEClust"), 
  lty = 'blank',
  cex = 2,
  cat.cex = 1.75,
  margin=0.2,
  fontfamily='sans',
  cat.fontfamily='sans')
  
dev.off()

# =======================================================================================
#test for differential expression of KRAUS1 and KRAUS3
DEseq<- nbinomTest(KRAUS_cds, "leaf1", "leaf3")

#rename columns and add info about species and comparison
DEseq["samples"]<-"13"
DEseq["species"]<-"KRAUS"
DEseq<-rename(DEseq, c("pval" = "DESeq_pval", "padj"="DESeq_padj"))

#read in results from DGEClust
DGEClust<-read.table(file="Bioinformatics/Selaginella/DGEClust/KRAUS13_pvals.txt", header=T)

#sort rows
DGEClust<-DGEClust[order(as.numeric(rownames(DGEClust))),]

#add DGESClust data to results
DEseq["DGEClust_pval"]<-DGEClust["pval"]
DEseq["DGEClust_padj"]<-DGEClust["padj"]

write.table(DEseq[DEseq$DESeq_padj < 0.1 & DEseq$DGEClust_padj < 0.1, ], 
  "Google Drive/Selaginella/DGEClust/KRAUS13_overlap.txt")

#remove missing values that will cause an error in the following step
DEseq<-na.omit(DEseq)

#this will replace extreme values with zero, which makes the plot nicer
DEseq[DEseq$DESeq_padj < 1e-5, ]["DESeq_padj"]<-0

#make scatterplot of adjusted p-values and Venn Diagram of overlap
pdf(file="Google Drive/Selaginella/DGEClust/KRAUS13.pdf")
ggplot(DEseq, aes(x=DESeq_padj, y=DGEClust_padj)) +
  geom_point(size=1.2) +
  scale_x_log10() +
  scale_y_log10() +
  geom_hline(yintercept=0.1, size=.2) +
  geom_vline(xintercept=0.1, size=.2) +
  theme_bw() +
  ggtitle("KRAUS13")
  
draw.pairwise.venn(area1=nrow(DEseq[DEseq$DGEClust_padj < 0.1, ]), 
  area2=nrow(DEseq[DEseq$DESeq_padj < 0.1, ]), 
  cross.area=nrow(DEseq[DEseq$DESeq_padj < 0.1 & DEseq$DGEClust_padj < 0.1, ]), 
  fill=c('blue', 'red'), 
  cat.col = c('blue', 'red'),
  category=c("DESeq", "DGEClust"), 
  lty = 'blank',
  cex = 2,
  cat.cex = 1.75,
  margin=0.2,
  fontfamily='sans',
  cat.fontfamily='sans')
  
dev.off()

# =======================================================================================
#test for differential expression of KRAUS1 and KRAUS4
DEseq<- nbinomTest(KRAUS_cds, "leaf1", "leaf4")

#rename columns and add info about species and comparison
DEseq["samples"]<-"14"
DEseq["species"]<-"KRAUS"
DEseq<-rename(DEseq, c("pval" = "DESeq_pval", "padj"="DESeq_padj"))

#read in results from DGEClust
DGEClust<-read.table(file="Bioinformatics/Selaginella/DGEClust/KRAUS14_pvals.txt", header=T)

#sort rows
DGEClust<-DGEClust[order(as.numeric(rownames(DGEClust))),]

#add DGESClust data to results
DEseq["DGEClust_pval"]<-DGEClust["pval"]
DEseq["DGEClust_padj"]<-DGEClust["padj"]

write.table(DEseq[DEseq$DESeq_padj < 0.1 & DEseq$DGEClust_padj < 0.1, ], 
  "Google Drive/Selaginella/DGEClust/KRAUS14_overlap.txt")

#remove missing values that will cause an error in the following step
DEseq<-na.omit(DEseq)

#this will replace extreme values with zero, which makes the plot nicer
DEseq[DEseq$DESeq_padj < 1e-5, ]["DESeq_padj"]<-0

#make scatterplot of adjusted p-values and Venn Diagram of overlap
pdf(file="Google Drive/Selaginella/DGEClust/KRAUS14.pdf")
ggplot(DEseq, aes(x=DESeq_padj, y=DGEClust_padj)) +
  geom_point(size=1.2) +
  scale_x_log10() +
  scale_y_log10() +
  geom_hline(yintercept=0.1, size=.2) +
  geom_vline(xintercept=0.1, size=.2) +
  theme_bw() +
  ggtitle("KRAUS14")
  
draw.pairwise.venn(area1=nrow(DEseq[DEseq$DGEClust_padj < 0.1, ]), 
  area2=nrow(DEseq[DEseq$DESeq_padj < 0.1, ]), 
  cross.area=nrow(DEseq[DEseq$DESeq_padj < 0.1 & DEseq$DGEClust_padj < 0.1, ]), 
  fill=c('blue', 'red'), 
  cat.col = c('blue', 'red'),
  category=c("DESeq", "DGEClust"), 
  lty = 'blank',
  cex = 2,
  cat.cex = 1.75,
  margin=0.2,
  fontfamily='sans',
  cat.fontfamily='sans')
  
dev.off()
# =======================================================================================
#test for differential expression of KRAUS2 and KRAUS3
DEseq<- nbinomTest(KRAUS_cds, "leaf2", "leaf3")

#rename columns and add info about species and comparison
DEseq["samples"]<-"23"
DEseq["species"]<-"KRAUS"
DEseq<-rename(DEseq, c("pval" = "DESeq_pval", "padj"="DESeq_padj"))

#read in results from DGEClust
DGEClust<-read.table(file="Bioinformatics/Selaginella/DGEClust/KRAUS23_pvals.txt", header=T)

#sort rows
DGEClust<-DGEClust[order(as.numeric(rownames(DGEClust))),]

#add DGESClust data to results
DEseq["DGEClust_pval"]<-DGEClust["pval"]
DEseq["DGEClust_padj"]<-DGEClust["padj"]

write.table(DEseq[DEseq$DESeq_padj < 0.1 & DEseq$DGEClust_padj < 0.1, ], 
  "Google Drive/Selaginella/DGEClust/KRAUS23_overlap.txt")

#remove missing values that will cause an error in the following step
DEseq<-na.omit(DEseq)

#this will replace extreme values with zero, which makes the plot nicer
DEseq[DEseq$DESeq_padj < 1e-5, ]["DESeq_padj"]<-0

#make scatterplot of adjusted p-values and Venn Diagram of overlap
pdf(file="Google Drive/Selaginella/DGEClust/KRAUS23.pdf")
ggplot(DEseq, aes(x=DESeq_padj, y=DGEClust_padj)) +
  geom_point(size=1.2) +
  scale_x_log10() +
  scale_y_log10() +
  geom_hline(yintercept=0.1, size=.2) +
  geom_vline(xintercept=0.1, size=.2) +
  theme_bw() +
  ggtitle("KRAUS23")
  
draw.pairwise.venn(area1=nrow(DEseq[DEseq$DGEClust_padj < 0.1, ]), 
  area2=nrow(DEseq[DEseq$DESeq_padj < 0.1, ]), 
  cross.area=nrow(DEseq[DEseq$DESeq_padj < 0.1 & DEseq$DGEClust_padj < 0.1, ]), 
  fill=c('blue', 'red'), 
  cat.col = c('blue', 'red'),
  category=c("DESeq", "DGEClust"), 
  lty = 'blank',
  cex = 2,
  cat.cex = 1.75,
  margin=0.2,
  fontfamily='sans',
  cat.fontfamily='sans')
  
dev.off()
