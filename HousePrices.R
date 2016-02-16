library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(scales)
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
                       
                       'V16'='', 
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
    print(url)
    Prices<-read.csv(url, header=F)
    AllPrices<-rbind(AllPrices, Prices)
}
names(AllPrices) <- c('unique_id',
                      'price_paid',
                      'deed_date',
                      'postcode',
                      'property_type',
                      'new_build',
                      'estate_type',
                      'saon',
                      'paon',
                      'street',
                      'locality',
                      'town',
                      'district',
                      'county',
                      'transaction_category',
                      'linked_data_uri'
                      )
AllPrices$street<- ifelse(AllPrices$paon == 10 & AllPrices$street == 'MARMADUKE STREET','10 MARMADUKE STREET', AllPrices$street)
AllPrices$Street<- ordered(AllPrices$street, levels=c('10 MARMADUKE STREET',
                                                       'MONTGOMERY STREET', 
                                                       'MARMADUKE STREET', 
                                                       'MONMOUTH STREET', 
                                                       'MERIONETH STREET', 
                                                       'MAIDSTONE STREET', 
                                                       'MARGATE STREET', 
                                                       'NOTTINGHAM STREET', 
                                                       'NEWPORT STREET'
                                                      )
                          )
filter(AllPrices, grepl('STREET', street)) %>% 
  ggplot(., aes(x=as.Date(deed_date, format="%Y-%m-%d"), y=as.numeric(price_paid))) +
  geom_point(aes(colour=Street)) +
  geom_smooth(method="loess", span=.3) +
  scale_colour_manual(values=c('red', rev(colorRampPalette(brewer.pal(9,"Blues"))(10)[3:10]))) +
  scale_x_date(breaks='year', labels = date_format("'%y")) +
  xlab("Date") +
  ylab("Price")
