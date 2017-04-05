
rm(list=ls())
setwd('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\expression_clustering')
library(data.table)
library(dplyr)

############################################
#Data process
############################################

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
exprinpath <- lapply(rlt,function(x)x[[1]])
geneinpath <- unique(unlist(lapply(rlt,function(x){colnames(x$data)})))
geneoutpath <- expr[,!colnames(expr)%in%geneinpath]

#PCA

library(GenABEL)
scale_ivn <- function(x){apply(x,2,rntransform)}
qpca <- function(A,rank=0,ifscale=F){
  if(!ifscale){A <- scale_ivn(A)}else{A <- scale(A)[,]}
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
rownames(x) <- colnames(expr)

############################################
#Distance
############################################

# save(x,file='x.rda')

############################################
#Summary
############################################

load('expression_clustering_label.rda')
rm(list=ls()[!ls()%in%ls(pattern='rlt.label|x|expr|scale_ivn|qpca')])

qpca2 <- function(A,ifscale){
  A.pca <- qpca(A,ifscale=ifscale)
  A.qpca <- qpca(A,rank=which(A.pca$prop>=0.9)[1],ifscale=ifscale)
  A.qpca$X[,1:which(A.qpca$prop>=0.9)[1],drop=F]
}
p_ginv_sq <- function(X,p){
  X.eigen = eigen(X);
  X.rank = sum(X.eigen$values>1e-8);
  X.value = X.eigen$values[1:X.rank]^(-1*p);
  if (length(X.value)==1){
    D = as.matrix(X.value);
  }else{
    D = diag(X.value);
  }
  rlt = X.eigen$vectors[,1:X.rank] %*% D %*% t(X.eigen$vectors[,1:X.rank]);
  return(rlt);
}
mrank <- function(X){
  X.svd = svd(X);
  X.rank = sum(X.svd$d>1e-6);
  return(X.rank);
}
mrank_sq <- function(X){
  X.eigen = eigen(X);
  X.rank = sum(Re(X.eigen$values)>1e-6);
  return(X.rank);
}
CCA_chisq_test <- function(rho,n,p,q){
  tstat = -1*n*sum(log(1-rho^2));
  p_value = pchisq(tstat,(p*q),lower.tail=FALSE);
  return(p_value);          
}
cca <- function(A,B){
  n = nrow(A);
  p = mrank(A);
  q = mrank(B);
  if (p <= q){
    X = A;
    Y = B;
  }else{
    X = B;
    Y = A;
  }
  R = p_ginv_sq(cov(X),0.5) %*% cov(X,Y) %*% p_ginv_sq(cov(Y),1) %*% cov(Y,X) %*% p_ginv_sq(cov(X),0.5);
  k = mrank_sq(R);
  d = Re(eigen(R)$values);
  rho = d[1:k]^(0.5);
  rho[rho >= 0.9999]=0.9;
  chisq_p = CCA_chisq_test(rho,n,p,q);
  return(c("chisq_p"=chisq_p,"df"=p*q));
}
ccap <- function(A,B){as.numeric(cca(A,B)[1])}

expr.pathout <- lapply(unique(rlt.label),function(i){
  x <- expr[,which(rlt.label==i),drop=F]
  qpca2(x,T)
})

expr2.pathout <- lapply(unique(rlt.label),function(i){
  x <- expr[,which(rlt.label==i),drop=F]
  x
})

i <- 0
rlt <- sapply(expr.pathout,function(x){
  print(i <<- i+1)
  sapply(expr.pathout,function(y){
    ccap(x,y)
  })
})
rlt2 <- rlt
# rlt2[rlt2 >= (0.1/length(rlt))] <- 1
rlt2 <- rlt*length(rlt)/2
rlt2 <- ifelse(rlt2>1,1,rlt2)

corrplot::corrplot(((1-rlt2)[,]))
heatmap(rlt2)

#Clustering Plot

sum(sapply(expr2.pathout,ncol))

test <- lapply(expr2.pathout,cor)
test <- lapply(test,colMeans)
test <- lapply(test,sort,decreasing=T)
test <- lapply(1:length(test),function(i){
  x <- expr2.pathout[[i]]
  sel <- test[[i]]
  sel <- names(sel)[1:min(length(sel),10)]
  x <- x[,colnames(x)%in%sel,drop=F]
  colnames(x) <- paste('group',i,colnames(x),sep="_")
  x
})

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
heatmap2 <- function(dat,col=rainbow(ncol(dat), start = 0, end = 1)){
  if(is.list(dat)){dat <- do.call(cbind,dat)}
  dat <- (as.matrix(dat))[,]
  # heatmap(t(dat),Colv=NA,col=col)
  # heatmap(t(dat),Colv=NA,Rowv=NA,col=col)
  heatmap(t(dat),Rowv=NA,col=col)
}

heatmap3 <- function(single_plot,p=T){
  if(is.list(single_plot)){single_plot <- do.call(cbind,single_plot)}
  if(p){single_plot <- 2*(pnorm(single_plot)-0.5)}
  heatmap(t(single_plot),Rowv=NA,
          col=colorRampPalette(c("green",'black',"red"))(1000),
          labCol=F,
          scale=c('row','column',"none"))
}
heatmap3(test[c(1:5,7,8,9)])
heatmap3(test[c(11,12,14,15,16,20,21)])
heatmap3(test[c(22:23,25,27,28,29,30)])

names(expr2.pathout) <- paste('cluster',1:length(expr2.pathout),sep="_")
expr_cluster <- c(exprinpath,expr2.pathout)
save(expr_cluster,file='expr_cluster.rda')

test <- lapply(expr_cluster[1:30],function(x) x[,1:min(10,ncol(x)),drop=F])
save(test,file='test.rda')
