
########################################
# Clustering Validation
########################################

rm(list=ls())
# devtools::install_github('wenrurumon/mysrc/lrm')
library(lrm)
setwd('C:\\Users\\zhu2\\Documents\\signaling\\codes\\')
source('sparse_2sem_final.R')
source('local_cnif_macro.R')
source('C:\\Users\\zhu2\\Documents\\signaling\\codes\\\\CNIF.R')
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
# Group Network
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
save(data.grp,rlt,file='C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\methylation_net\\methylation_network_rlt.rda')

########################################
# inGroup Network
########################################

rm(list=ls())
setwd('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\methylation_net')
source('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\source_grpsem.R')
source('C:\\Users\\zhu2\\Documents\\signaling\\codes\\\\CNIF.R')
load('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\methylation_net\\methylation_network_rlt.rda')
plotnet <- function(x,mode='directed',main=NULL){
  diag(x) <- 0
  plot(graph_from_adjacency_matrix(t(x),mode=mode),
       edge.arrow.size=.3,
       vertex.size=3,
       vertex.label.cex=1,
       edge.width=.1,
       main=main)
}

rlt2.model <- function(i,j,door=2){
  datai <- data.grp[[i]][[j]]
  semi <- sparse_2sem(Y=datai,lambda=0.5)
  gc()
  if(door==2) {
    dagi <- try(CNIF(data=datai,init.adj=semi[[1]],max_parent = 3))
  } else {
    dagi <- NA
  }
  if(!is.matrix(dagi)){
    gc()
    dagi <- try(CNIF(data=datai,init.adj=semi[[1]],max_parent = 2))
  }
  try(plotnet(dagi))
  return(dagi)
}

rlt2.model2 <- function(i,j,door=2){
  datai <- data.grp[[i]][[j]]
  test <- qpca(datai)
  test <- test$X[,1:which(test$prop>=0.3)[1],drop=F]
  datai <- lm(datai~test)$residuals
  semi <- sparse_2sem(Y=datai,lambda=0.25)
  gc()
  if(door==2) {
    dagi <- try(CNIF(data=datai,init.adj=semi[[1]],max_parent = 3))
  } else {
    dagi <- NA
  }
  if(!is.matrix(dagi)){
    gc()
    dagi <- try(CNIF(data=datai,init.adj=semi[[1]],max_parent = 2))
  }
  try(plotnet(dagi,main=paste(i,j,door)))
  return(dagi)
}
# rlt2.model2(1,1)

# rlt2 <- lapply(1:2,function(i){
#   lapply(1:length(data.grp[[i]]),function(j){
#     rlt <- try(rlt2.model(i,j))
#     rlt
# })})

setwd('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\methylation_net\\')
# load('methylation_site_network.rda')
load('methylation_site_network_new.rda')
length(rlt2)
# rlt2 <- list()
k <- 0
for(i in (length(rlt2)+1):22){
  print(i)
  rlt <- lapply(1:length(data.grp[[i]]),function(j){
    k <<- j
    print('\n');print(paste('model i=',i,'j=',j));print('\n')
    if(j%in%(c(11,35)+1)){
      rlt <- try(rlt2.model2(i,j,1))
    } else {
      rlt <- try(rlt2.model2(i,j))  
    }
    rlt
  })
  rlt2[[i]] <- rlt
  save(rlt2,file='methylation_site_network_new.rda')
}
names(rlt2) <- names(data.grp)
save(rlt2,data.grp,file='methylation_site_network_new.rda')
