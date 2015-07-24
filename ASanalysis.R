# Heath O'Brien (heath.obrien-at-gmail-dot-com)
# AS analysis
# 16 July 2015

# last modified 22 July 2015

# Plotting and analysis of survey data

library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(grid)
require(ordinal) #clm

setwd("~/Documents/R")
ASdata <- read.delim("ASdata.txt")
summary(ASdata)
#makes a comment
#next line: replace w/
ASdata$Gender<-ifelse(ASdata$Q1 == 2, "Male", ifelse(ASdata$Q1 == 1, "Female", "NA"))

#go ahead and remove responses from people who's gender is NA
ASdata<-subset(ASdata, Gender != 'NA')
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

ASdata$Q16<-ifelse(ASdata$Q16 < 3, "Agree", 
                  ifelse(ASdata$Q16 == 3, "Slightly agree", 
                         ifelse(ASdata$Q16 == 4, "Neutral",
                                ifelse(ASdata$Q16 == 5, "Slightly disagree", 
                                       ifelse(ASdata$Q16 > 5, "Disagree", "NA")))))
ASdata$Q16<-factor(ASdata$Q16, levels=c("Agree", "Slightly agree", "Neutral", "Slightly disagree", "Disagree"))
ggplot(subset(ASdata, Q16 != 'NA'), aes(Q16, fill=Gender)) + 
  geom_bar(position='dodge') +
  facet_grid(Position ~ .) +
  xlab("response") +
  ggtitle('I feel I have a good work-life balance') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = unit(c(1,1,1.1,1), "cm")) +
  scale_fill_manual(values=(brewer.pal(12, "Paired")[c(1,2)])) +
  scale_x_discrete(labels=c("Agree", "Slightly agree", "Neutral", "Slightly disagree", "Disagree"))
grid.text("© Heath O'Brien 2015", x=unit(.99, "npc"), y=unit(.01, "npc"), just=c("right", "bottom"))
lapply(ASdata[, c("Q16", "Position", "Gender")], table)
ftable(xtabs(~ Q16 + Position + Gender, data = ASdata))

MOD.1 <- clm(Q16 ~ Position + Gender, link="logit", data = ASdata)
summary(MOD.1)
count <- ASdata %>% 
           select(Gender, Position, Q16) %>% 
           group_by(Gender, Position, Q16) %>% 
           tally()
count <- ASdata %>% 
           select(Position, Gender, Q16) %>%           
           filter(complete.cases(.)) %>% 
           expand(Position, Gender, Q16) %>% 
           left_join(count)
count[is.na(count)] <- 0
count <- count %>% 
           group_by(Gender, Position) %>% 
           summarise(Total = sum(n)) %>% left_join(count)

newdat = data.frame(Position = c("Academic", "Staff", "Student", "Academic", "Staff"), 
                    Gender=c("Male", "Female", "Male", "Female", "Male"))
count<-cbind(count, predict(MOD.1, count, type = "prob", se.fit=TRUE))

ggplot(count, aes(Q16, n/Total, fill=Gender)) + 
  geom_bar(stat="identity", position='dodge') +
  facet_grid(Position ~ .) +
  xlab("response") +
  ggtitle('I feel I have a good work-life balance') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = unit(c(1,1,1.1,1), "cm")) +
  scale_fill_manual(values=(brewer.pal(12, "Paired")[c(1,2)])) +
  geom_point(aes(Q16, fit))
grid.text("© Heath O'Brien 2015", x=unit(.99, "npc"), y=unit(.01, "npc"), just=c("right", "bottom"))
