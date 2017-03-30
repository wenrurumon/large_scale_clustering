
rm(list=ls())
setwd("/home/zhu/rushdata/expression_clustering")
load('gene_feature.rda')

cmax <- cor(t(x)); gc();
dmax <- as.matrix(dist(x)); gc();

library(igraph)
g.cmax <- igraph::graph_from_adjacency_matrix(cmax>0.7,mode='undirected')
g.dmax <- igraph::graph_from_adjacency_matrix(dmax<0.16,mode='undirected')

system.time(rw.cmax <- cluster_walktrap(g.cmax))
system.time(rw.dmax <- cluster_walktrap(g.dmax))
