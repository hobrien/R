# Heath O'Brien (heath.obrien-at-gmail-dot-com)
# AS analysis
# 16 July 2015

# last modified 22 July 2015

# Plotting and analysis of survey data

library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(grid)

setwd("~/Documents/R")
ASdata <- read.delim("ASdata.txt")
summary(ASdata)
#makes a comment
#next line: replace w/
ASdata$Gender<-ifelse(ASdata$Q1 == 2, "Male", ifelse(ASdata$Q1 == 1, "Female", "NA"))
ASdata$Position<-ifelse(ASdata$Q2 <= 3, "Academic", ifelse(ASdata$Q2 == 6, "Student", "Staff"))
#Replace answer with label in any question, here Q7
ASdata$Q7<-ifelse(ASdata$Q7 == 1, "7.3", 
                  ifelse(ASdata$Q7 == 2, "7.3-10", 
                  ifelse(ASdata$Q7 == 3, "10-12",
                  ifelse(ASdata$Q7 == 4, "12+", "NA"))))
#stacked bar plot of seniority versus gender
ggplot(ASdata, aes(factor(Q2), fill=Gender))+geom_bar()
#plots hours worked versus gender and seniority
ggplot(subset(ASdata, Gender != 'NA'), aes(factor(Q9), fill=Gender)) + 
  geom_bar() +
  facet_grid(Position ~ Gender) +
  xlab("response") +
  ggtitle('In order to do my job effectively, I need to work this many hours per day') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = unit(c(1,1,1.1,1), "cm")) +
  scale_fill_manual(values=(brewer.pal(12, "Paired")[c(1,2)]))
grid.text("© Heath O'Brien 2015", x=unit(.99, "npc"), y=unit(.01, "npc"), just=c("right", "bottom"))
