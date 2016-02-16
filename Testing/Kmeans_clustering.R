library(ggplot2)
library(reshape2)

#Make Parallel coordinates plot of movie ratings
#This will give a plot of the percentage of votes for each rating (from 1 to 10) for the most rated 840 movies.

#This is from Hadley's Use R book 

popular <- subset(movies, votes > 1e4)
ratings <- popular[, 7:16]
ratings$.row <- rownames(ratings)
molten <- melt(ratings, id = ".row")
ggplot(molten, aes(variable, value, group = .row)) +
  geom_line(alpha = 1 / 20, position=position_jitter(w=0, h=2))

#Add k-means clustering of movies with similar ratings
k = 12
cl <- kmeans(ratings[1:10], k)
ratings$cluster <- reorder(factor(cl$cluster), popular$rating)
levels(ratings$cluster) <- seq_along(levels(ratings$cluster)) #I'm pretty unclear on what this does
molten <- melt(ratings, id = c(".row", "cluster"))
ggplot(molten, aes(variable, value, group = .row, colour=cluster)) +
  geom_line(alpha = 1 / 20, position=position_jitter(w=0, h=2))

#can add faceting to see the groupings easier
ggplot(molten, aes(variable, value, group = .row, colour=cluster)) +
  geom_line(alpha = 1 / 20, position=position_jitter(w=0, h=2)) +
  facet_wrap(~ cluster)
