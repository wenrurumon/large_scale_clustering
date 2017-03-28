
rm(list=ls())
setwd('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\expression_clustering')
library(data.table)
library(dplyr)

#Data process

dir()
map <- fread("RNA_Seq_gene_expression_levels_info.txt")
raw <- fread("RNA_Seq_gene_expression_levels_clean.txt")
pid <- raw$ProjectID
raw <- data.matrix(dplyr::select(raw,-Sample,-ProjectID))
colnames(raw) <- map$genenm
rownames(raw) <- pid
gene <- colnames(raw)
raw[pid=='82317494',] <- rbind(colMeans(raw[pid=='82317494',]),colMeans(raw[pid=='82317494',]))
raw <- raw[-which(pid=='82317494')[2],]
raw <- sapply(unique(gene),function(g){
  rowMeans(raw[,gene==g,drop=F])
})
expr <- scale(raw[,(apply(raw,2,var)>0),drop=F])[,]

load('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\summary\\expr_network_20170308.rda')
geneinpath <- unique(unlist(lapply(rlt,function(x){colnames(x$data)})))
geneoutpath <- expr[,!colnames(expr)%in%geneinpath]

#PCA

library(GenABEL)
scale_ivn <- function(x){apply(x,2,rntransform)}
qpca <- function(A,rank=0){
  A <- scale_ivn(A)
  A.svd <- svd(A)
  if(rank==0){
    d <- A.svd$d
  } else {
    d <- A.svd$d-A.svd$d[min(rank+1,nrow(A),ncol(A))]
  }
  d <- d[d > 1e-10]
  r <- length(d)
  prop <- d^2; prop <- cumsum(prop/sum(prop))
  d <- diag(d,length(d),length(d))
  u <- A.svd$u[,1:r,drop=F]
  v <- A.svd$v[,1:r,drop=F]
  x <- u%*%sqrt(d)
  y <- sqrt(d)%*%t(v)
  z <- x %*% y
  rlt <- list(rank=r,X=x,Y=y,Z=x%*%y,prop=prop)
  return(rlt)
}
A <- t(expr)
system.time(A.pca <- qpca(A))
system.time(A.qpca <- qpca(A,rank=which(A.pca$prop>=0.9)[1]))
x <- A.qpca$X[,1:which(A.qpca$prop>=0.9)[1],drop=F]
system.time(xdist <- dist(x))
save(xdist,file='xdist.rda')
