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
  if (nrow(DEseq[DEseq$DESeq_padj < 1e-5, ]["DESeq_padj"]) > 0 ) {
    DEseq[DEseq$DESeq_padj < 1e-5, ]["DESeq_padj"] <- 0
  }
  return(DEseq)
}

#read in raw count data
selaginellaCountTable =  read.table( file="Bioinformatics/Selaginella/Counts/Selag_counts.txt" , header=TRUE, row.names=1 )

for (species in c('KRAUS', 'UNC', 'WILD', 'MOEL')) {
  if (species == 'MOEL') {
     max=3 
     condition = factor( c( "leaf1", "leaf2", "leaf3" ) )
     CountTable = selaginellaCountTable [c(paste(species, 1, sep=''), paste(species, 2, sep=''), paste(species, 3, sep=''))]
  }  else {
     max=4
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
       textfile<-paste("Google Drive/Selaginella/DGEClust/", name, "_overlap.txt", sep="")
       plotfile<-paste("Google Drive/Selaginella/DGEClust/", name, "_overlap.pdf", sep="")  
       results<-merge_results(DESeq.results, DGEClust, name, textfile)
       pdf(plotfile)
       print(ggplot(results, aes(x=DESeq_padj, y=DGEClust_padj)) +
          geom_point(size=1.2) +
          scale_x_log10() +
          scale_y_log10() +
          geom_hline(yintercept=0.1, size=.2) +
          geom_vline(xintercept=0.1, size=.2) +
          theme_bw() +
          ggtitle(name))
       if ( nrow(results[results$DGEClust_padj < 0.1, ]) > 0 && nrow(results[results$DESeq_padj < 0.1, ]) > 0 ) {
         draw.pairwise.venn(area1=nrow(results[results$DGEClust_padj < 0.1, ]), 
            area2=nrow(results[results$DESeq_padj < 0.1, ]), 
            cross.area=nrow(results[results$DESeq_padj < 0.1 & results$DGEClust_padj < 0.1, ]), 
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
       dev.off()
    }
  }
}
