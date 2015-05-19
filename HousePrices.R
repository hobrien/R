library(ggplot2)
library(dplyr)
library(RColorBrewer)
AllPrices <-data.frame('V1'='', 
                       'V2'='', 
                       'V3'='',
                       'V4'='', 
                       'V5'='', 
                       'V6'='', 
                       'V7'='', 
                       'V8'='', 
                       'V9'='',
                       'V10'='', 
                       'V11'='',
                       'V12'='',
                       'V13'='', 
                       'V14'='', 
                       'V15'='', 
                       stringsAsFactors=FALSE
                       ) 

for (Street in c('Montgomery', 
                'Marmaduke', 
                'Monmouth', 
                'Merioneth', 
                'Maidstone', 
                'Margate', 
                'Nottingham', 
                'Newport'
                )) {
  print(Street)
    #to look up all houses: paste("http://landregistry.data.gov.uk/app/ppd/?et[]=lrcommon%3Afreehold&et[]=lrcommon%3Aleasehold&limit=100&max_date=", format(Sys.Date(), "%d+%B+%Y"), "&min_date=1+January+1999&nb[]=true&nb[]=false&ptype[]=lrcommon%3Adetached&ptype[]=lrcommon%3Asemi-detached&ptype[]=lrcommon%3Aterraced&ptype[]=lrcommon%3Aflat-maisonette&street=", Street, "&town=Bristol", sep="")
    url <- paste("http://landregistry.data.gov.uk/app/ppd/ppd_data.csv?et%5B%5D=lrcommon%3Afreehold",
                 "limit=all",
                 "min_date=1+January+1999",
                 paste("max_date", format(Sys.Date(), "%d+%B+%Y"), sep='='),
                 "nb%5B%5D=false",
                 "ptype%5B%5D=lrcommon%3Aterraced",
                 paste("street", Street, sep='='),
                 "town=Bristol", 
                 sep="&"
                )
    Prices<-read.csv(url, header=F)
    AllPrices<-rbind(AllPrices, Prices)
}
names(AllPrices) <- c('ID',
                      'Price', 
                      'Date', 
                      'Post Code', 
                      'Terrace', 
                      'New', 
                      'Leasehold',
                      'Blank1',
                      'Number',
                      'Street',
                      'City1',
                      'City2',
                      'City3',
                      'City4',
                      'URL'
                      )
AllPrices$Street<- ordered(AllPrices$Street, levels=c('MONTGOMERY STREET', 
                                                       'MARMADUKE STREET', 
                                                       'MONMOUTH STREET', 
                                                       'MERIONETH STREET', 
                                                       'MAIDSTONE STREET', 
                                                       'MARGATE STREET', 
                                                       'NOTTINGHAM STREET', 
                                                       'NEWPORT STREET'
                                                      )
                          )
colorRampPalette(brewer.pal(9,"Blues"))(100)
filter(AllPrices, grepl('STREET', Street)) %>% 
  ggplot(., aes(x=as.Date(Date, format="%Y-%m-%d"), y=as.numeric(Price))) +
  geom_point(aes(colour=Street)) +
  geom_smooth(method="loess") +
  scale_colour_manual(values=rev(colorRampPalette(brewer.pal(9,"Blues"))(10)[3:10])) +
  xlab("Date") +
  ylab("Price")
