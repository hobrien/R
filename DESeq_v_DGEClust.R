library(ggplot2)
library(plyr)
library(DESeq)
library(VennDiagram)

setwd("/Users/HeathOBrien/")

merge_results<-function(DEseq, DGEclust, name, textfile){
  #rename columns and add info about species and comparison
  DEseq["samples"]<-name
  DEseq<-rename(DEseq, c("pval" = "DESeq_pval", "padj"="DESeq_padj"))

  #sort rows
  DGEClust<-DGEClust[order(as.numeric(rownames(DGEClust))),]

  #add DGESClust data to results
  DEseq["DGEClust_pval"]<-DGEClust["pval"]
  DEseq["DGEClust_padj"]<-DGEClust["padj"]

  #write list of overlaps to textfile
  write.table(DEseq[DEseq$DESeq_padj < 0.1 & DEseq$DGEClust_padj < 0.1, ], textfile)

  #remove missing values that will cause an error in the following step
  DEseq<-na.omit(DEseq)

  #this will replace extreme values with zero, which makes the plot nicer
  DEseq[DEseq$DESeq_padj < 1e-5, ]["DESeq_padj"]<-0
  return(DEseq)
}
scatterplot<-function(DEseq, name) {
  #make scatterplot of adjusted p-values and Venn Diagram of overlap
  sp =  ggplot(DEseq, aes(x=DESeq_padj, y=DGEClust_padj)) +
    geom_point(size=1.2) +
    scale_x_log10() +
    scale_y_log10() +
    geom_hline(yintercept=0.1, size=.2) +
    geom_vline(xintercept=0.1, size=.2) +
    theme_bw() +
    ggtitle(name)
    return(sp)
}
vennplot <-function(DEseq) {
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
 }


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
#I'm not sure how best to loop this
KRAUS12 <- nbinomTest(KRAUS_cds, "leaf1", "leaf2")
#KRAUS13 <- nbinomTest(KRAUS_cds, "leaf1", "leaf3")
#KRAUS14 <- nbinomTest(KRAUS_cds, "leaf1", "leaf4")
#KRAUS23 <- nbinomTest(KRAUS_cds, "leaf2", "leaf3")
#KRAUS24 <- nbinomTest(KRAUS_cds, "leaf2", "leaf4")
#KRAUS34 <- nbinomTest(KRAUS_cds, "leaf3", "leaf4")

DGEClust<-read.table(file="Bioinformatics/Selaginella/DGEClust/KRAUS12_pvals.txt", header=T)
textfile<-"Google Drive/Selaginella/DGEClust/KRAUS12_overlap.txt"
plotfile<-"Google Drive/Selaginella/DGEClust/KRAUS12_overlap.pdf"  
name = "KRAUS12"

KRAUS12<-merge_results(KRAUS12, DGEClust, name, textfile)
pdf(plotfile)
scatterplot(KRAUS12, name)
vennplot(KRAUS12)
dev.off()
