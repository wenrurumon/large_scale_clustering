rm(list=ls())

library(ggplot2)
library(ggdendro)

hc <- hclust(dist(USArrests), "ave")
ggdendrogram(hc, rotate = FALSE, size = 2)
hc$labels[hc$order]

dat <- scale(as.matrix(USArrests))

hc <- heatmap(dat)
dat.hc <- dat[hc$rowInd,hc$colInd]
dat.m <- melt(dat,id=rownames(dat))
colnames(dat.m) <- c('city','case','value')
dat.m <- ddply(dat.m,.(case),mutate,score=scale(value))

ggplot(data = dat.m, aes(x = case, y = city))+
  geom_tile(aes(fill = score),colour='white') +
  scale_fill_gradient(low = "white",high = "steelblue")
