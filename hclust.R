
rm(list=ls())

#Server

load('xdist.rda')
cluster <- function(xdist,members=NULL){
  if(is.null){
   xdist2 <- xdist  
  }
  xdist2 <- xdist[rownames(xdist)%in%members,colnames(xdist)%in%members,drop=F]
  xhc <- hclust(xdist2)
  return(xhc)
}

test1 <- cluster(xdist)
