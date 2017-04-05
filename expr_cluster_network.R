
setwd('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\expression_clustering')
load('expr_cluster.rda')
load('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\disease.rda')
library(data.table)
library(dplyr)
# 
# library(GenABEL)
# scale_ivn <- function(x){apply(x,2,rntransform)}
# qpca <- function(A,rank=0,ifscale=F){
#   if(!ifscale){A <- scale_ivn(A)}else{A <- scale(A)[,]}
#   A.svd <- svd(A)
#   if(rank==0){
#     d <- A.svd$d
#   } else {
#     d <- A.svd$d-A.svd$d[min(rank+1,nrow(A),ncol(A))]
#   }
#   d <- d[d > 1e-10]
#   r <- length(d)
#   prop <- d^2; prop <- cumsum(prop/sum(prop))
#   d <- diag(d,length(d),length(d))
#   u <- A.svd$u[,1:r,drop=F]
#   v <- A.svd$v[,1:r,drop=F]
#   x <- u%*%sqrt(d)
#   y <- sqrt(d)%*%t(v)
#   z <- x %*% y
#   rlt <- list(rank=r,X=x,Y=y,Z=x%*%y,prop=prop)
#   return(rlt)
# }
# 
# qpca2 <- function(A,ifscale){
#   if(ncol(A)==1){return(cbind(scale(A)[,]))}
#   A.pca <- qpca(A,ifscale=ifscale)
#   A.qpca <- qpca(A,rank=which(A.pca$prop>=0.9)[1],ifscale=ifscale)
#   A.qpca$X[,1:which(A.qpca$prop>=0.9)[1],drop=F]
# }
# 
# expr_cluster <- lapply(expr_cluster,function(x){
#   x <- x[match(disease$projid,rownames(x)),,drop=F]
#   x[,apply(x,2,var)>0,drop=F]
# })
# 
# expr_cluster <- lapply(expr_cluster,function(x){
#   qpca2(x,ifscale=T)
# })

# i <- 0
# ccap_mat <- sapply(expr_cluster2,function(x){
#   print(i<<-i+1)
#   sapply(expr_cluster2,function(y){
#     ccap(x,y)
#   })
# })

#######################
#SEM
#######################

source('C:\\Users\\zhu2\\Documents\\sample_grouplasso_network\\scr\\grpsem.R')
system.time(test <- grpsem(expr_cluster[1:5],lambda=0.7))
plot(graph_from_adjacency_matrix(t(test)))

########################

library(flare)
Y <- expr_cluster

groups <- seq(1,662,20)
groups[length(groups)] <- 662
groups <- lapply(2:length(groups)-1,function(i){groups[[i]]:groups[[i+1]]})

#grpslim
tryout <- function(sel){
  Ymat.datai <- data4sem(y=y,x=x[sel])
  test <- slim(X=Ymat.datai$x,Y=Ymat.datai$
                 y,lambda=0.5)
  tapply(test$beta0[[1]]!=0,Ymat.datai$index,any)
}

test.slim <- lapply(1:5,function(i){
  print(i)
  y <- Y[[i]]
  x <- Y[-i]
  out <- lapply(groups[1:5],function(sel){
    try(tryout(sel))
  })
  
  

})
test.slim <- do.call(rbind,test.slim)+0

#########################################

setwd('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\expression_clustering')
load('expr_cluster.rda')
load('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\disease.rda')
library(data.table)
library(dplyr)

expr_cluster2 <- sapply(expr_cluster,rowMeans)
dim(expr_cluster2)

test <- sapply(1:663,function(i){
  print(i)
  test <- flare::slim(X=expr_cluster2[,-i],Y=expr_cluster2[,i],lambda=0.2)
  out <- rep(FALSE,663)
  out[-i] <- (test$beta!=0)
  out
})
plotnet <- function(x,mode='undirected'){
  diag(x) <- 0
  plot(graph_from_adjacency_matrix(t(x),mode=mode),
       edge.arrow.size=.1,
       vertex.size=3,
       vertex.label.cex=1,
       edge.width=.1)
}
plotnet(test)
undnet633 <- test
save(undnet633,file='undnet633.rda')
