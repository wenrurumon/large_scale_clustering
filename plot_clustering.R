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

##################

rm(list=ls())
dat <- USArrests

ggheatmap <- function(dat,dimnames=NULL,col=c('white','white','steelblue')){
  dat <- scale(as.matrix(dat))[,]
  dat.m <- melt(dat,id=rownames(dat))
  if(!is.null(dimnames)){
    colnames(dat.m)[1:2] <- dimnames
  } else {
    colnames(dat.m)[1:2] <- c('id','variables')
  }
  s <- paste0('ggplot(data=dat.m,aes(x=',colnames(dat.m)[2],',y=',colnames(dat.m)[1],'))+',
         'geom_tile(aes(fill=value),colour="',col[1],'")+',
         'scale_fill_gradient(low="',col[2],'",high="',col[3],'")')
  eval(parse(text=s))
}

ggheatmap(dat,col=c('gray','white','red'))

heatmap(t(single_plot),Rowv=NA,
        col=colorRampPalette(c("green",'black',"red"))(80),
        labCol=F,
        scale=c('row','column',"none"),
        margins=c(7,4))

