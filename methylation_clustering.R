
rm(list=ls())
setwd('/home/zhu/rushdata/processed/methylation')
qpca <- function(A,rank=0,ifscale=TRUE){
  if(ifscale){A <- scale(as.matrix(A))[,]}
  A.svd <- svd(A)
  if(rank==0){
    d <- A.svd$d
  } else {
    d <- A.svd$d-A.svd$d[min(rank+1,nrow(A),ncol(A))]
  }
  d <- d[d > 1e-8]
  r <- length(d)
  prop <- d^2; info <- sum(prop)/sum(A.svd$d^2);prop <- cumsum(prop/sum(prop))
  d <- diag(d,length(d),length(d))
  u <- A.svd$u[,1:r,drop=F]
  v <- A.svd$v[,1:r,drop=F]
  x <- u%*%sqrt(d)
  y <- sqrt(d)%*%t(v)
  z <- x %*% y
  rlt <- list(rank=r,X=x,Y=y,Z=x%*%y,prop=prop,info=info)
  return(rlt)
}
qpca2 <- function (A, ifscale = T, l1 = 0.9, l2 = 0.9) 
{
    A.pca <- qpca(A)
    A.qpca <- qpca(A, rank = which(A.pca$prop >= l1)[1])
    A.qpca$X[, 1:which(A.qpca$prop >= l2)[1], drop = F]
}
f <- dir()
data <- lapply(f,function(x){
  print(x)
  load(x)
  raw <- scale(sapply(pmethy,function(x){
   rowMeans(x$methy)
  }))[,]
  raw2 <- t(qpca2(t(raw)))
  colnames(raw2) <- colnames(raw)
  cmax <- cor(raw2)
  dmax <- as.matrix(dist(t(raw2)))
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
save(rlt,file='methylation_clustering.rda')



