
#################################
# Project
#################################

rm(list=ls())

setwd('C:\\Users\\zhu2\\Documents\\signaling\\codes\\')
source('sparse_2sem_final.R')
source('local_cnif_macro.R')
source('CNIF.R')
sourceCpp("score_function_regression.cpp")
sourceCpp("simple_cycle.cpp")
sourceCpp("initial_sem.cpp")
setwd('C:\\Users\\zhu2\\Documents\\sample_grouplasso_network')
sourceCpp("scr\\score_function_regression.cpp")
sourceCpp("scr\\simple_cycle.cpp")
sourceCpp("scr\\initial_sem.cpp")
source("scr\\local_cnif_macro.R")
source("scr\\grpsem.R")
library(flare)
library(grplasso)
library(data.table)
library(dplyr)
library(igraph)
load("C:/Users/zhu2/Documents/getpathway/model20170215/expression_clustering/expr_cluster_list.rda")
load('C:/Users/zhu2/Documents/getpathway/model20170215/expression_clustering/expnet_outpath.rda')
load('C:/Users/zhu2/Documents/getpathway/model20170215/rlt_p2pfinal.rda')

qpca2 <- function(x){
  x <- x[,apply(x,2,var)>0,drop=F]
  if(ncol(x)==1){return(x)}
  x1 <- which(qpca(x)$prop>=0.9)[1]
  if(x1<ncol(x)){
    x2 <- qpca(x,rank=x1)
  } else {
    x2 <- qpca(x)
  }
  x2$X[,1:which(x2$prop>=0.9)[1],drop=F]
}
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
plotnet <- function(x,mode='undirected'){
  diag(x) <- 0
  plot(graph_from_adjacency_matrix(t(x),mode=mode),
       edge.arrow.size=.1,
       vertex.size=3,
       vertex.label.cex=1,
       edge.width=.1)
}

#################################
# Connect all clusters
#################################

allpaths <- lapply(do.call(c,exprincluster),function(x){
  x[match(rownames(exprincluster[[1]][[1]]),rownames(x)),,drop=F]
})
# i <- 0
allpaths.qpca <- sapply(allpaths,function(x){
  # print(i<<-i+1)
  qpca2(x)[,1]
})
system.time(cluster.network <- sparse_2sem(Y=allpaths.qpca,lambda=0.4,times=10))
plotnet((cluster.network[[1]]>=0.8)[-1:-299,-1:-299])
# colnames(cluster.network[[1]][1:299,1:299])
# colnames(rlt_p2pfinal[[2]])[1:299]

cluster.network[[1]][match(colnames(rlt_p2pfinal[[2]]),colnames(cluster.network[[1]])),match(colnames(rlt_p2pfinal[[2]]),colnames(cluster.network[[1]]))] <- rlt_p2pfinal[[2]]
# plotnet(cluster.network[[1]][1:299,1:299])
# plotnet(rlt_p2pfinal[[2]])
rlt_p2pfinal[[2]][1,rlt_p2pfinal[[2]][1,,drop=F]==1,drop=F]
# which(colnames(cluster.network[[1]])=='Metabolic pathways')
cluster.network[[1]][299,which(cluster.network[[1]][299,,drop=F]>0),drop=F]
cluster_network <- (cluster.network[[1]]>=.8)+0
cluster_CNIF <- CNIF(data=allpaths.qpca,init.adj=cluster_network,max_parent=2)
cluster_CNIF[match(colnames(rlt_p2pfinal[[2]]),colnames(cluster_CNIF)),match(colnames(rlt_p2pfinal[[2]]),colnames(cluster_CNIF))] <- rlt_p2pfinal[[2]]
is.dag(igraph::graph_from_adjacency_matrix(cluster_CNIF))
sum(colSums(cluster_CNIF)>0|rowSums(cluster_CNIF)>0)

##########################
#Prepare Output
##########################

exprincluster <- (allpaths)
cluster_network <- cluster_CNIF
load("C:/Users/zhu2/Documents/getpathway/model20170215/expression_clustering/temp2.rda")
load("C:/Users/zhu2/Documents/getpathway/model20170215/summary/expnet.rda")
expnet <- c(lapply(rlt2,function(x)x[[1]]),rlt)
dim2 <- function(x){
  x1 <- nrow(x)
  x2 <- ncol(x)
  x1 <- ifelse(is.null(x1),NA,x1)
  x2 <- ifelse(is.null(x2),NA,x2)
  c(x1,x2)
}
test <- cbind(t(sapply(allpaths,dim2)),t(sapply(expnet,dim2)))
save(exprincluster,cluster_network,expnet,file='C:/Users/zhu2/Documents/getpathway/model20170215/rlt_cluster_network.rda')

###########################
#cluster 2 cluster
###########################

rlt_all_split <- lapply(1:length(allpaths),function(i){
  print(i)
  Yi <- allpaths[[i]]
  Xi.list <- allpaths[cluster_network[i,]>0]
  Yinet <- expnet[[i]][,1:nrow(expnet[[i]]),drop=F]
  if(length(Xi.list)==0){
    semi <- Yinet
  } else {
    semi <- sparse_2sem(Y=Yi,Y.fixed=Yinet,X=do.call(cbind,Xi.list),lambda=0.5)[[1]]  
  }
  semi
})
save(exprincluster,cluster_network,rlt_all_split,expnet,file='C:/Users/zhu2/Documents/getpathway/model20170215/rlt_cluster_network.rda')

