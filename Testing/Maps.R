#see https://github.com/hadley/ggplot2/wiki/plotting-polygon-shapefiles
#ward map from http://www.arcgis.com/home/item.html?id=cd5721b36d39413f82f60913a1f4823d
library(maptools)
library(ggplot2)
library(rgdal)
library(plyr)
library(rgeos)  #gSimplify

#url <- "http://www.nyc.gov/html/dcp/download/bytes/nybb_15b.zip"
#fil <- basename(url)
#if (!file.exists(fil)) download.file(url, fil)
#fils <- unzip(fil)

gor <- readShapeSpatial('BristolWards/BristolWards.shp')
gor@data$id = rownames(gor@data)
#can simplify map with gor <- gSimplify(gor, 0.05), but it isn't necessary in this case
gor.points <- fortify(gor, region = "id")
gor.df=join(gor.points, gor@data, by="id")
ggplot(gor.df) + 
  aes(long,lat,group=group,fill=WardName) + 
  geom_polygon() +
  geom_path(color="white") +
  coord_equal() +
  #scale_fill_brewer("Electoral Wards") +
  labs(x="",y="") +
  scale_x_continuous(breaks=NULL) +
  scale_y_continuous(breaks=NULL) 


