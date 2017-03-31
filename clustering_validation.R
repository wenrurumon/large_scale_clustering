###############
#summary

rm(list=ls())
setwd('/home/zhu/rushdata/expression_clustering/')
load('methylation_cmax_dmax.rda')
load('methylation_clustering.rda')
i <- 1

xdist <- data[[i]]$dmax
xcluster <- rlt[[i]]$xtree
xcluster <- do.call(rbind,lapply(1:length(xcluster),function(i){
  cbind(i,xcluster[[i]])
}))
xcluster <- as.numeric(xcluster[match(colnames(xdist),xcluster[,2]),1])
