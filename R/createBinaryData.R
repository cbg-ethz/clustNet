createBinaryData <- function(clusterSizes, n, seed=101) {

# generate binary data for each cluster and combine the matrices

  set.seed(seed)

  kclust<-length(clusterSizes)

  for(ii in 1:kclust){

    oneclustsize<-clusterSizes[ii]

    binarycluster<-matrix(0,nrow=oneclustsize,ncol=n)
    binarytest<-matrix(0,nrow=oneclustsize,ncol=n)

    for(kk in 1:n){ # build data for this cluster
      beta<-runif(1)
      binarycluster[,kk]<-sample.int(2,oneclustsize,replace=TRUE,prob=c(beta,1-beta))-1
    }

    if(ii==1){ # combine the clusters
      binaryclusterdata<-binarycluster
    } else {
      binaryclusterdata<-rbind(binaryclusterdata,binarycluster)
    }

  }

  return(binaryclusterdata)
}
