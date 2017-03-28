
rm(list=ls())

#Server

load('xdist.rda')
#xhc <- hclust(xdist)

cluster <- function(xdist,members=NULL){
  if(is.matrix(xdist)){
    mdist <- xdist
    xdist <- as.dist(xdist)
  } else {
    mdist <- as.matrix(xdist)
  }
  if(is.null(members)){
    xdist2 <- xdist  
  } else {
    xdist2 <- as.dist(mdist[attr(xdist,"Labels")%in%members,attr(xdist,"Labels")%in%members,drop=F])
  }
  xhc <- hclust(xdist2)
  xtree <- cutree(xhc,2)
  xtree <- lapply(1:2,function(i){
    names(which(xtree==i))
   })
  xdists <- lapply(xtree,function(x){
    mdist[rownames(mdist)%in%x,colnames(mdist)%in%x,drop=F]
  })
  return(list(xtree=xtree,xdists=xdists))
}

system.time(rlt <- cluster(xdist))
i <- which(sapply(rlt$xtree,length)==max(sapply(rlt$xtree,length)))
system.time(rltsub <- cluster(rlt$xdists[[i]],members=rlt$xtree[[i]]))
