#hacked together from http://www.uni-kiel.de/psychologie/rexrepos/posts/regressionOrdinal.html
#and http://www.ats.ucla.edu/stat/r/dae/ologit.htm

#Ordinal Logit Regression
require(foreign) #read.dta
require(ggplot2)
require(MASS) #polr
require(ordinal) #clm

dat <- read.dta("http://www.ats.ucla.edu/stat/data/ologit.dta")
head(dat)
lapply(dat[, c("apply", "pared", "public")], table)
ftable(xtabs(~ public + apply + pared, data = dat))

ggplot(dat, aes(x = apply, y = gpa)) +
  geom_boxplot(size = .75) +
  geom_jitter(alpha = .5) +
  facet_grid(pared ~ public, margins = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

#I can't seem to get se.fit from MASS
#m <- polr(apply ~ pared + public + gpa, data = dat, Hess=TRUE)

#clm works a trick tho
clmFit <- clm(apply ~ pared + public + gpa, link="logit", data = dat)
summary(m)
newdat <- data.frame(
  pared = rep(0:1, 200),
  public = rep(0:1, each = 200),
  gpa = rep(seq(from = 1.9, to = 4, length.out = 100), 4))

newdat <- cbind(newdat, predict(clmFit, newdat, type = "probs", se.fit=TRUE))
