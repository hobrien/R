library(ggplot2)
library(plyr)
library(DESeq)
library(VennDiagram)

setwd("/Users/HeathOBrien/")

MergeResults<-function(DEseq, DGEclust, name, textfile){
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
  if (nrow(DEseq[DEseq$DESeq_padj < 1e-5, ]["DESeq_padj"]) > 0 ) {
    DEseq[DEseq$DESeq_padj < 1e-5, ]["DESeq_padj"] <- 0
  }
  return(DEseq)
}
scatterplot<-function(pvalues, DESeq_cutoff, DGEClust_cutoff, name) {
  #make scatterplot of adjusted p-values and Venn Diagram of overlap
  sp = ggplot(pvalues, aes(x=DESeq_padj, y=DGEClust_padj)) +
    geom_point(size=1.2) +
    scale_x_log10() +
    scale_y_log10() +
    geom_hline(yintercept=DGEClust_cutoff, size=.2) +
    geom_vline(xintercept=DESeq_cutoff, size=.2) +
    theme_bw() +
    ggtitle(name)
    return(sp)
}
vennplot <-function(pvalues, DESeq_cutoff, DGEClust_cutoff) {
   if ( nrow(pvalues[pvalues$DGEClust_padj < DGEClust_cutoff, ]) > 0 && nrow(pvalues[pvalues$DESeq_padj < DESeq_cutoff, ]) > 0 ) {
      if ( nrow(pvalues[pvalues$DGEClust_padj < DGEClust_cutoff, ]) < nrow(pvalues[pvalues$DESeq_padj < DESeq_cutoff, ]) ) {
         draw.pairwise.venn(area1=nrow(pvalues[results$DGEClust_padj < DGEClust_cutoff, ]), 
            area2=nrow(pvalues[pvalues$DESeq_padj < DESeq_cutoff, ]), 
            cross.area=nrow(pvalues[pvalues$DESeq_padj < DESeq_cutoff & pvalues$DGEClust_padj < DGEClust_cutoff, ]), 
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
      else {
         draw.pairwise.venn(area1=nrow(pvalues[pvalues$DESeq_padj < DESeq_cutoff, ]), 
            area2=nrow(pvalues[pvalues$DGEClust_padj < DGEClust_cutoff, ]), 
            cross.area=nrow(pvalues[pvalues$DESeq_padj < DESeq_cutoff & pvalues$DGEClust_padj < DGEClust_cutoff, ]), 
            fill=c('red', 'blue'), 
            cat.col = c('red', 'blue'),
            category=c("DGEClust", "DESeq"), 
            lty = 'blank',
            cex = 2,
            cat.cex = 1.75,
            margin=0.2,
            fontfamily='sans',
            cat.fontfamily='sans')
      }
   }
}
 
 
#read in raw count data
selaginellaCountTable =  read.table( file="Bioinformatics/Selaginella/Counts/Selag_counts.txt" , header=TRUE, row.names=1 )
species.list <- c('KRAUS', 'MOEL', 'UNC', 'WILD')
DESeq_cutoff = 0.1
for (DGEClust_cutoff in c(0.1, 0.01, 0.001)) {
   plotfile<-paste("Google Drive/Selaginella/DGEClust/Overlap_", DGEClust_cutoff, ".pdf", sep="")  
   pdf(plotfile)
   #Analyse Species comparisons
   condition = factor( c( 'KRAUS', 'KRAUS', 'KRAUS', 'KRAUS', 'MOEL', 'MOEL', 'MOEL', 'UNC', 'UNC', 'UNC', 'UNC', 'WILD', 'WILD', 'WILD', 'WILD' ) )
   DESeq<-newCountDataSet( selaginellaCountTable, condition)
   DESeq <- estimateSizeFactors(DESeq)
   DESeq <- estimateDispersions(DESeq)
   for (x in 1:3) {
      for (y in (x+1):4) {
         sample1 <- species.list[x]
         sample2 <-species.list[y]
         name = paste(sample1, sample2, sep="_")
         print(name)
         DESeq.results <- nbinomTest(DESeq, sample1, sample2)
         DGEClust<-read.table(file=paste("Bioinformatics/Selaginella/DGEClust/pvals_", name, ".txt", sep=""), header=T)
         DGEClust<-rename(DGEClust, c("Posteriors" = "pval", "FDR" = "padj"))
         textfile<-paste("Google Drive/Selaginella/DGEClust/", name, "_overlap_", DGEClust_cutoff, ".txt", sep="")
         results<-MergeResults(DESeq.results, DGEClust, name, textfile)
         print(scatterplot(results, DESeq_cutoff, DGEClust_cutoff, name))
         print(vennplot(results, DESeq_cutoff, DGEClust_cutoff))
      }
   }
   for (species in species.list) {
      if (species == 'MOEL') {
         max=3 
         condition = factor( c( "leaf1", "leaf2", "leaf3" ) )
         CountTable = selaginellaCountTable [c(paste(species, 1, sep=''), paste(species, 2, sep=''), paste(species, 3, sep=''))]
      }
      else {
         max=2#4
         condition = factor( c( "leaf1", "leaf2", "leaf3", "leaf4" ) )
         CountTable = selaginellaCountTable[c(paste(species, 1, sep=''), paste(species, 2, sep=''), paste(species, 3, sep=''), paste(species, 4, sep=''))]
      }

      DESeq<-newCountDataSet( CountTable, condition)
      DESeq <- estimateSizeFactors(DESeq)
      DESeq <- estimateDispersions(DESeq, method="blind", sharingMode="fit-only")
      for (sample1 in 1:(max-1)){
         for (sample2 in (sample1+1):max){
            DESeq.results <- nbinomTest(DESeq, paste("leaf", sample1, sep=""), paste("leaf", sample2, sep=""))
            name = paste(species, sample1, sample2, sep="")
            DGEClust<-read.table(file=paste("Bioinformatics/Selaginella/DGEClust/", name, "_pvals.txt", sep=""), header=T)
            textfile<-paste("Google Drive/Selaginella/DGEClust/", name, "_overlap_", DGEClust_cutoff, ".txt", sep="")
            print(scatterplot(results, DESeq_cutoff, DGEClust_cutoff, name))
            print(vennplot(results, DESeq_cutoff, DGEClust_cutoff))
         }
      }
   }
   dev.off()
}       
    
