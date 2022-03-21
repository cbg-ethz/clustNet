checkmembership <- function(kclust,kclusttrue,truelabels,estmemb) {
  relabelmatrix<-matrix(nrow=kclust,ncol=kclusttrue) #i row_num label, j col_num estmemb
  estlabels<-list()
  for (i in 1:kclust) {
    estlabels[[i]]<-which(estmemb==i)
  }
  for (i in 1:kclust) {
    for (j in 1:kclusttrue) {
      relabelmatrix[i,j]<-length(which(estlabels[[i]]%in%truelabels[[j]]))
    }
  }
  rowcol<-clue::solve_LSAP(relabelmatrix,TRUE)
  res<-list()
  res$relabel<-as.vector(rowcol)
  res$ncorr<-0
  for (j in 1:min(kclust,kclusttrue)) {
    res$ncorr<-res$ncorr+relabelmatrix[j,res$relabel[j]]
  }
  res$relabelmatrix<-relabelmatrix
  return(res)
}

#' @importFrom stats runif
generatetriple<-function(n) {
  resmat<-matrix(nrow=n,ncol=3)
  for (i in 1:n) {
    ordr<-sample.int(3,3)
    res<-vector(length=3)
    res[1]<-runif(1,min=0,max=1)
    res[2]<-runif(1,min=0,max=1-res[1])
    res[3]<-1-res[1]-res[2]
    resmat[i,ordr]<-res
  }
  return(resmat)
}

#' @importFrom stats runif
generatefour<-function(n) {
  resmat<-matrix(nrow=n,ncol=4)
  for (i in 1:n) {
    ordr<-sample.int(4,4)
    res<-vector(length=4)
    res[1]<-runif(1,min=0,max=1)
    res[2]<-runif(1,min=0,max=1-res[1])
    res[3]<-runif(1,min=0,max=1-res[1]-res[2])
    res[4]<-1-res[1]-res[2]-res[3]
    resmat[i,ordr]<-res
  }
  return(resmat)
}

# propersample <- function(x){if(length(x)==1) x else sample(x,1)}

# calcloglike <- function(samplescores,tau) {
#  # samplescores<-t(samplescores)
#   maxscorey<-apply(samplescores,1,max) # find the max of each column
#   loglike<-sum(log(colSums(t(exp(samplescores-maxscorey))*tau))+maxscorey) # remove max for numerical stability and exponentiate
#   return(loglike)
# }

# reassignsamples <- function(samplescores,numsamps){
#   newclustermembership <-rep(0,numsamps) # to store the new cluster
#   for(s in 1:numsamps){ # run through the samples
#     clusterscores<-samplescores[s,]
#     maxscorey<-max(clusterscores) # take the maximum
#     maxscoreelem<-which(clusterscores==maxscorey)
#     newclustermembership[s]<-propersample(maxscoreelem)
#   }
#   return(newclustermembership)
# }

# relativeprobs <- function(samplescores,numsamps){
#   relativeprobabs <-rep(0,numsamps) # to store the relative probabilities
#   for(s in 1:numsamps){ # run through the samples
#     clusterscores<-samplescores[s,]
#     maxscorey<-max(clusterscores) # take the maximum
#     shifty<-exp(clusterscores-maxscorey)
#     rescaley<-shifty/sum(shifty)
#     relativeprobabs[s]<-max(rescaley) # relative probabilities
#   }
#   return(relativeprobabs)
# }

# allrelativeprobs <- function(samplescores,numsamps){
#   relativeprobabs <-samplescores # to store the relative probabilities
#   for(s in 1:numsamps){ # run through the samples
#     clusterscores<-samplescores[s,]
#     maxscorey<-max(clusterscores) # take the maximum
#     shifty<-exp(clusterscores-maxscorey)
#     rescaley<-shifty/sum(shifty)
#     relativeprobabs[s,]<-rescaley # relative probabilities
#   }
#   return(relativeprobabs)
# }

relativeprobswithtau <- function(sampleprobs,tau){
  temp<-tau*t(sampleprobs)
  relativeprobabswithtau<-1/colSums(temp)*t(temp) # ugly code
  return(relativeprobabswithtau)
}

avescore <- function(samplescores,numsamps){
  averagescores <-rep(0,numsamps) # to store the relative probabilities
  for(s in 1:numsamps){ # run through the samples
    clusterscores<-samplescores[s,]
    maxscorey<-max(clusterscores) # take the maximum
    shifty<-exp(clusterscores-maxscorey) # exponentiate
    rescaley<-log(mean(shifty))+maxscorey # mean and turn back to log
    averagescores[s]<-rescaley # averagescore
  }
  return(averagescores)
}
reassignsamplesprop <- function(samplescores,numsamps,gamma){
  newclustermembership <-rep(0,numsamps) # to store the new cluster
  for(s in 1:numsamps){ # run through the samples
    clusterscores<-samplescores[s,]*gamma
    maxscorey<-max(clusterscores) # take the maximum
    shifty<-exp(clusterscores-maxscorey)
    rescaley<-shifty/sum(shifty)
    scorelength<-length(rescaley)
    newclustermembership[s]<-sample.int(scorelength,1,prob=rescaley) # sample according to scores
  }
  return(newclustermembership)
}


#' @title netCluster
#'
#' @description Network-based clustering
#'
#' @param myData Data to be clustered
#' @param kclust Number of clusters
#' @param nbg Number of covariates
#' @param itLim Maximum number of iterations
#' @param EMseeds Seeds
#' @param BBMMClust Binary clustering before network-based clustering (TRUE by default)
#' @param edgepmat Matrix of penalized edges in the search space
#' @param bdepar Hyperparameters for structure learning (BDE score)
#'
#' @return a list containing the clusterMemberships and "assignprogress"
#' @export
netCluster <- function(myData,kclust=3,nbg=0,itLim=20, EMseeds=1, BBMMClust=TRUE, edgepmat=NULL, bdepar=list(chi = 0.5, edgepf = 16)){

  # Binary clustering
  startseed <- EMseeds[1]
  nIterations <- 10
  chi <- 0.5 # pseudocounts for the Beta prior

  if(BBMMClust){
    binClust <- BBMMclusterEM(binaryMatrix = myData, chi = chi, kclust = kclust, startseed = startseed, nIterations = nIterations, verbose=TRUE)
    newallrelativeprobabs <- binClust$relativeweights
    newclustermembership <- binClust$newclustermembership
    newclustermembership <- reassignsamples(newallrelativeprobabs)
  }

  # # number of iterations (and respective seeds)
  # nSeeds <- 2
  # # number of covariate variables
  # nbg <- 3
  ss <- dim(myData)[1]
  # #number of clusters
  # kclust<-3
  # # iteration limit
  # itLim <- 20

  # total number of variables
  nn <- dim(myData)[2]
  #number of variables (without covariates)
  n<-nn-nbg

  # #prior pseudo counts
  # chixi<-0.5
  #define when EM converges
  err<-1e-6
  # #edge penalization factor
  # edgepfy<-16
  #define different seeds to run EM
  # EMseeds<-c(101,102,103,104,105)
  # EMseeds<-c(100)+c(1:nSeeds)
  #number of EM attempts to ensure highest likelihood
  EMn<-length(EMseeds)
  #number of iterations in the internal cycle
  nit.internal<-10


  #to store accuracy of cluster centers
  # centerprogress<-list()
  #accuracy of assignments for different EM runs
  assignprogress<-list()
  #to store scores against clusters for each EM run
  scoresprogress<-list()
  #to store cluster probabilities for all sample for each EM run
  probs<-list()
  #to store relabelling
  relabs<-list()
  #to store cluster centers
  clustercenters<-list()
  #to store scores against clusters
  scoresagainstclusters<-matrix(ncol=kclust,nrow=ss)

  for (s in 1:EMn) {
    diffy<-1
    cnt<-1
    assignprogress[[s]]<-list()
    assignprogress_local<-list()
    # assignprogress_local$corrvec<-numeric()
    assignprogress_local$likel<-numeric()
    set.seed(EMseeds[s])
    print(paste("EM seed",EMseeds[s]))
    # centers<-list()
    for (i in 1:kclust) {
      # centers[[i]]<-as.data.frame(matrix(ncol=3))
      clustercenters[[i]]<-matrix(rep(0,(n+nbg)^2),nrow=n+nbg)
    }


    if(!BBMMClust){
      #generate random assignment of belonging to each cluster for each sample
      newallrelativeprobabs<-generatetriple(ss)
      newclustermembership<-reassignsamples(newallrelativeprobabs)
    }


    # #learn how many samples were asssigned correctly
    # res<-checkmembership(kclust,kclusttrue,truelabels,newclustermembership)
    # print(paste("number of correctly assigned samples by random assignment:", res$ncorr,
    #             "of total", ss, "samples"))

    # EM cycle
    while (diffy>err&cnt<itLim) {
      allrelativeprobabs<-newallrelativeprobabs
      allprobprev<-newallrelativeprobabs
      coltots<-colSums(allrelativeprobabs) + chixi # add prior to clustersizes
      tauvec<-coltots/sum(coltots)

      #outer EM cycle, learning cluster centers
      #define cluster centers and assign probabilities
      for (k in 1:kclust) {

        # #learn background nodes
        # for (i_nbg in 1:nbg){
        #   for (i_n in 1:n){
        #     #chi squared test
        #     testRes <- wtd.chi.sq(myData[,i_n],myData[,n+i_nbg],weight=allrelativeprobabs[,k])
        #     siglevel <- 0.05/(kclust+n)
        #     #add to initial DAG structure
        #     if (testRes[3]<siglevel){
        #       clustercenters[[k]][n+i_nbg,i_n] <- 1
        #     }else{
        #       clustercenters[[k]][n+i_nbg,i_n] <- 0
        #     }
        #   }
        # }

        #define score parameters
        if (nbg>0){
          scorepar<-BiDAG::scoreparameters("bde",myData, edgepmat = edgepmat,
                                    weightvector=allrelativeprobabs[,k],
                                    bdepar=bdepar, bgnodes=(n+1):(n+nbg))
        }else{
          scorepar<-BiDAG::scoreparameters("bde",myData, edgepmat = edgepmat,
                                    weightvector=allrelativeprobabs[,k],
                                    bdepar=bdepar)
        }

        #find MAP DAG using iterative order MCMC
        maxfit<-BiDAG::iterativeMCMC(scorepar,addspace=clustercenters[[k]],verbose=FALSE)
        #store it
        clustercenters[[k]]<-maxfit$DAG

      }

      #internal EM cycle, estimating parameters
      for (i in 1:nit.internal) {
        allrelativeprobabs<-newallrelativeprobabs
        coltots<-colSums(allrelativeprobabs) + chixi # add prior to clustersizes
        tauvec<-coltots/sum(coltots)
        for (k in 1:kclust) {
          if (nbg>0){
            scorepar<-BiDAG::scoreparameters("bde",as.data.frame(myData), edgepmat = edgepmat,
                                             weightvector=allrelativeprobabs[,k],
                                             bdepar=bdepar, bgnodes=(n+1):(n+nbg))
          }else{
            scorepar<-BiDAG::scoreparameters("bde",as.data.frame(myData), edgepmat = edgepmat,
                                             weightvector=allrelativeprobabs[,k],
                                             bdepar=bdepar
          }
          scorepar$n <- n # to avoid to scoring over background nodes
          scoresagainstclusters[,k]<-BiDAG::scoreagainstDAG(scorepar,clustercenters[[k]])
          scorepar$n <- n+nbg # recet after scoring
        }
        newallrelativeprobabsnotau<-allrelativeprobs(scoresagainstclusters)
        newallrelativeprobabs<-relativeprobswithtau(newallrelativeprobabsnotau,tauvec)
      }

      diffy<-sum((allprobprev-newallrelativeprobabs)^2)
      newclustermembership<-reassignsamples(newallrelativeprobabs)
      # res<-checkmembership(kclust,kclusttrue=3,truelabels,newclustermembership)
      # assignprogress_local$corrvec[cnt]<-res$ncorr
      assignprogress_local$likel[cnt]<-calcloglike(scoresagainstclusters,tauvec)
      # for (k in 1:3) {
      #   centers[[k]][cnt,]<-compareDAGs(m2graph(clustercenters[[k]]),
      #                                   BNsBG[[res$relabel[k]]]$DAG)[c("SHD","TP","FP")]
      # }
      cnt<-cnt+1
    }
    print(paste("EM converged after",cnt-1,"iterations"))
    # print(paste("number of correctly assigned samples:",
    # res$ncorr, "of", ss))
    #store number of correct assignments and likelihood
    assignprogress[[s]]<-assignprogress_local
    #store center progress
    # centerprogress[[s]]<-centers
    scoresprogress[[s]]<-scoresagainstclusters
    probs[[s]]<-newallrelativeprobabs
    # relabs[[s]]<-res$relabel
  }

  return(list("clustermembership"=newclustermembership,"assignprogress"=assignprogress, "DAGs"=clustercenters))
}


getBestSeed <- function(assignprogress){
  #get the vector with loglikelihood reached at the EM convergence
  likelihoodvector<-unlist(lapply(assignprogress,function(x)x$likel[length(x$likel)]))
  # assignvector<-unlist(lapply(assignprogress,function(x)x$corrvec[length(x$likel)]))

  EMseeds <- 1:length(assignprogress)
  EMn<-length(EMseeds)

  #order loglikelihood vector
  so<-order(likelihoodvector,decreasing = TRUE)
  # summary.res<-cbind(likelihoodvector[so],round(assignvector[so]/ss,3),EMseeds[so])
  # colnames(summary.res)<-c("logL","acc","seed")
  # print(summary.res)

  #define best seed and corresponding accuracy
  # perccorr<-round(assignvector[so][1]/ss*100)
  bestseed<-EMseeds[so][1]
  # percworst<-round(assignvector[so][EMn]/ss*100)
  worstseed<-EMseeds[so][EMn]

  return(bestseed)
}


#' @title netClusterParallel
#'
#' @description Network-based clustering of multiple seeds using parallel computing
#'
#' @param myData Data to be clustered
#' @param kclust Number of clusters
#' @param nbg Number of covariates
#' @param itLim Maximum number of iterations
#' @param EMseeds seeds
#' @param BBMMClust binary clustering before network-based clustering (TRUE by default)
#' @param edgepmat a matrix of penalized edges in the search space
#' @param bdepar Hyperparameters for structure learning (BDE score)
#'
#' @return a list containing the clusterMemberships and "assignprogress"
#' @export
#'
netClusterParallel <- function(myData,kclust=3,nbg=0,itLim=20, EMseeds=1:5, BBMMClust=TRUE, edgepmat=NULL, bdepar=list(chi = 0.5, edgepf = 16)){

  # parallel computing of clustering
  nSeeds <- length(EMseeds)
  clusterResAll <- parallel::mclapply(1:nSeeds, function(i) {
    print(paste("Clustering iteration", i, "of", nSeeds))
    clusterRes <- netCluster(myData=myData,kclust=kclust,nbg=nbg,itLim=itLim, EMseeds=EMseeds[i], BBMMClust=BBMMClust, edgepmat=edgepmat, bdepar=bdepar)
    return(clusterRes)
  }, mc.cores = nSeeds)

  # get best performing seed
  bestSeed <- getBestSeed(lapply(clusterResAll, function(x) x[[2]][[1]]))
  bestRes <- clusterResAll[[bestSeed]]

  # return results of best seed
  return(list("clustermembership"=bestRes$clustermembership,"assignprogress"=bestRes$assignprogress,"DAGs"=bestRes$DAGs,"bestSeed"=bestSeed))
}

