
rm(list=ls())
setwd('/home/zhu/rushdata/processed/methylation')
f <- dir()
data <- lapply(f,function(x){
  print(x)
  load(x)
  raw <- scale(sapply(pmethy,function(x){
   rowMeans(x$methy)
  }))[,]
  cmax <- cor(raw)
  dmax <- as.matrix(dist(t(raw)))
  list(raw=raw,cmax=cmax,dmax=dmax)                    
})

setwd('/home/zhu/rushdata/expression_clustering')
save(data,file='methylation_cmax_dmax.rda')

#########################################

rm(list=ls())
load('/home/zhu/rushdata/expression_clustering/methylation_cmax_dmax.rda')
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

cluster2 <- function(xdist){
  if(is.matrix(xdist)){
    rlt <- cluster(xdist)
  }else{
    rlt <- xdist
  }  
  while(max(sapply(rlt$xtree,length))>200){
    i <- which(sapply(rlt$xtree,length)==max(sapply(rlt$xtree,length)))[1]
    system.time(rltsub <- cluster(rlt$xdists[[i]],members=rlt$xtree[[i]]))
    rlt$xtree <- c(rlt$xtree[-i],rltsub$xtree)
    rlt$xdists <- c(rlt$xdists[-i],rltsub$xdists)
    print(sapply(rlt$xtree,length))
  }
  return(rlt)
 }

###############

system.time(test <- cluster2(1-data[[1]]$cmax))
i <- 0
rlt <- lapply(data,function(x){
  print(i<<-i+1)
  return(cluster2(1-x$cmax))
})
