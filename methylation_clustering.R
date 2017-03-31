
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
