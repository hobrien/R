library("ggplot2")
library("ape")
library("plyr")
library("reshape")

CalculateDistances<-function(file){
   aln<-read.dna(file=file, format="fasta")
   dist<-dist.dna(aln, model="raw", pairwise.deletion=T)  #caluculate pairwise distances
   dist<-as.matrix(dist)  #convert to matrix
   dist[upper.tri(dist, diag=T)]<-NA  #convert values on and above diagonal to NAs to avoid double counting
   dist<-melt(dist)  #convert to dataframe
   dist<-dist[!is.na(dist$value),]   #removing NAs
   dist$Comparison=paste(sub("((comp)|_).*", "", dist$X1), sub("((comp)|_).*", "", dist$X2), sep="_")  #concatinate spcecies names (removing gene info)
   dist<-adply(dist, 1, transform, Comparison = paste(sort(strsplit(Comparison, "_")[[1]]), collapse="_"))  #sort species in comparison
   dist<-rename(dist, c("value" = "Distance"))
   dist$Locus<-sub("^([^.]*).*", "\\1", basename(file))  #add locus info from file name
   return(dist)
}

#setwd("/Users/HeathOBrien/Bioinformatics/Selaginella/SingleCopy/")
x<- data.frame(Comparison=c(), Distance=c(), Locus=c())

for (file in list.files('.')) {
   dist <- tryCatch(CalculateDistances(file), error = function(e) NULL)
   x<- rbind(x, dist[c("Distance", "Comparison", "Locus")])
}
#setwd("/Users/HeathOBrien/Google Drive/Selaginella/SingleCopy/")
pdf("Histograms.pdf")
ggplot(x, aes(Distance))+geom_histogram()+facet_wrap(~Comparison)+theme_bw()
dev.off()