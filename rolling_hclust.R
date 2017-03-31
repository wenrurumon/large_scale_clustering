
rm(list=ls())

#Server

load('x.rda')
system.time(xdist <- 1-cor(t(x)))
#load('xdist.rda')
#load('temp.rda')
#xhc <- hclust(xdist)

cluster <- function(xdist,members=NULL,speed=2){
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
  xtree <- cutree(xhc,speed)
  xtree <- lapply(1:speed,function(i){
    names(which(xtree==i))
   })
  xdists <- lapply(xtree,function(x){
    mdist[rownames(mdist)%in%x,colnames(mdist)%in%x,drop=F]
  })
  return(list(xtree=xtree,xdists=xdists))
}

rlt <- cluster(xdist)
while(max(sapply(rlt$xtree,length))>200){
  i <- which(sapply(rlt$xtree,length)==max(sapply(rlt$xtree,length)))[1]
  system.time(rltsub <- cluster(rlt$xdists[[i]],members=rlt$xtree[[i]]))
  rlt$xtree <- c(rlt$xtree[-i],rltsub$xtree)
  rlt$xdists <- c(rlt$xdists[-i],rltsub$xdists)
  print(sapply(rlt$xtree,length))
}

