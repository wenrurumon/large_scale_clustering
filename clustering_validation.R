
rm(list=ls())
setwd('/home/zhu/rushdata/expression_clustering/')
load('methylation_cmax_dmax.rda')
load('methylation_clustering.rda')

cluster_label <- function(xdist,xcluster){
  xcluster <- do.call(rbind,lapply(1:length(xcluster),function(i){
    cbind(i,xcluster[[i]])
  }))
  xcluster <- as.numeric(xcluster[match(colnames(xdist),xcluster[,2]),1])
  xcluster
}

cluster_validation <- function(xdist,xcluster){
  sapply(unique(xcluster),function(j){
    js <- which(xcluster==j)
    c(mean(xdist[js,js]),mean(xdist[js,-js]))
  })
}

# validation the clustering
validation <- lapply(1:length(data),function(i){
  print(i)
  xcluster <- cluster_label(data[[i]]$dmax,rlt[[i]]$xtree)
  rlt.vali_d <- cluster_validation(data[[i]]$dmax,xcluster)
  rlt.vali_c <- cluster_validation(data[[i]]$cmax,xcluster)
  list(xcluster,rlt.vali_d,rlt.vali_c)
})

save(validation,file='methylation_clustering_validation.rda')
###############
#export
rm(list=ls())

setwd('/home/zhu/rushdata/expression_clustering')
load('methylation_clustering.rda')
load('methylation_clustering_validation.rda')
load('methylation_cmax_dmax.rda')
load('pid.rda')

data.grp <- lapply(1:22,function(i){
	print(i)
	grpi <- validation[[i]][[1]]
	datai <- data[[i]]$raw
	datai <- datai[match(pid,rownames(datai)),,drop=F]
	datai <- datai[,apply(datai,2,var)>0,drop=F]
	colnames(datai) <- paste(
		gsub('[^0-9]','',names(data))[i],
		colnames(datai),
		sep=":"
	)
	data.grp <- lapply(unique(grpi),function(j){
			datai[,grpi==j,drop=F]
		})
	names(data.grp) <- paste0(
			gsub('[^0-9]','',names(data))[i],
			':cluster',1:length(data.grp)
		)
	return(data.grp)
	})
names(data.grp) <- names(data)
save(data.grp,file='methylation_clusters_data.rda')
