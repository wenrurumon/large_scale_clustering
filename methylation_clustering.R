
rm(list=ls())
setwd('/home/zhu/rushdata/processed/methylation')
f <- dir()
data <- lapply(f,function(x){
  load(x)
  sapply(pmethy,function(x){
   rowMeans(x$methy)
  })
})
