#Server

load('xdist.rda')
cluster <- function(xdist,members){
  xdist2 <- xdist[rownames(xdist)%in%members,colnames(xdist)%in%members,drop=F]
  xhc <- hclust(xdist2)
  return(xhc)
}
