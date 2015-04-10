library(ggplot2)
library(lubridate)
library(scales)

y<-read.table("~/Desktop/Desktop_files/transactions.txt", sep="\t")

y[y$V4> -800 & y$V4<800,]$V3<-NA  #get rid of labels for small transactions

y$V1 = as.Date(y$V1, "%d-%b-%y") #convert first column to dates

#calculate savings per month since lowest point

#this assumes that the low point for the current month occurred within the last 20 days
thisMonth <- seq(Sys.Date(), length = 2, by = "-20 days")[2]

CurrentLow = y[y$V5 == min(y[y$V1 > thisMonth, ]$V5),]$V5
CurrentLowDate = y[y$V5 == min(y[y$V1 > thisMonth, ]$V5),]$V1
MaxLow = y[y$V5 == min(y[y$V1 > "2013-06-01", ]$V5),]$V5
MaxLowDate = y[y$V5 == min(y[y$V1 > "2013-06-01", ]$V5),]$V1
NumMonths = month(CurrentLowDate)-month(MaxLowDate)+24  # Need to add 12 if start is before current year (this will have to be increased to 24 next month)

require(zoo)
pdf(paste("~/Desktop/Desktop_files/transactions", format(as.yearmon(Sys.Date()), "%b%y"), '.pdf', sep=''), width = 16, height = 8)
ggplot(y, aes(x=V1, y=V5, group=1, label=V3)) +
       geom_line(cex=.05) +
       geom_segment(aes(x=MaxLowDate, y=MaxLow, xend=CurrentLowDate, yend=CurrentLow), colour="royalblue", cex=.1) +
       geom_smooth(se=FALSE, colour="royalblue")+geom_point() +
       scale_colour_manual(values=c("red", "red", "black", "red", "red", "red", "red", "red")) +
       scale_y_continuous("Balance", breaks=0:10*1000)+scale_x_date("Date", breaks = "month", labels = date_format("%m-%y")) +
       geom_text(cex=3, hjust=0, aes(colour=V2))+theme(legend.position="none") +
       annotate("text", x=CurrentLowDate, y=250, label=sprintf("£%.2f / month", (CurrentLow-MaxLow)/NumMonths), hjust=1)
dev.off()