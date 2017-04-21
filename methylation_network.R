
rm(list=ls())

# devtools::install_github('wenrurumon/mysrc/lrm')
library(lrm)
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

setwd('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\expression_clustering\\')
load('methylation_clusters_data.rda')

datai <- do.call(cbind,lapply(data.grp[[15]][1:5],function(x){
  x <- x[,1:20]
  x
}))
single_plot <- datai
test <- heatmap(t(pnorm(scale(single_plot))),
        Rowv=NA,
        # Colv=NA,
        col=colorRampPalette(c("green",'black',"red"))(80),
        labCol=F,
        scale=c('row','column',"none"),
        margins=c(7,4))

########################################

rm(list=ls())
source('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\source_grpsem.R')
load('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\expression_clustering\\methylation_clusters_data.rda')
plotnet <- function(x,mode='directed'){
  diag(x) <- 0
  plot(graph_from_adjacency_matrix(t(x),mode=mode),
       edge.arrow.size=.3,
       vertex.size=3,
       vertex.label.cex=1,
       edge.width=.1)
}

i <- 1
datai <- data.grp[[i]]

# x <- lapply(datai,qpca2,l1=.8,l2=.99)
# x <- lapply(x,function(x){
#   list(score=x)
# })
test <- model(datai,prop=.64,lambda=0.7,max.parent=3)
plotnet(test)
is.dag(graph_from_adjacency_matrix(t(test)))

rlt <- lapply(1:22,function(i){
  print(names(data.grp)[i])
  datai <- data.grp[[i]]
  gc();
  test <- try(model(datai,prop=.64,lambda=0.7,max.parent=4))
  if(!is.matrix(test)){
    test <- try(model(datai,prop=.64,lambda=0.7,max.parent=3))
  }
  plotnet(test)
  test
})

for(i in 1:22){plotnet(rlt[[i]])}

########################################
