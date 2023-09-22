checkmembership <- function(k_clust,kclusttrue,truelabels,estmemb) {
  relabelmatrix<-matrix(nrow=k_clust,ncol=kclusttrue) #i row_num label, j col_num estmemb
  estlabels<-list()
  for (i in 1:k_clust) {
    estlabels[[i]]<-which(estmemb==i)
  }
  for (i in 1:k_clust) {
    for (j in 1:kclusttrue) {
      relabelmatrix[i,j]<-length(which(estlabels[[i]]%in%truelabels[[j]]))
    }
  }
  rowcol<-clue::solve_LSAP(relabelmatrix,TRUE)
  res<-list()
  res$relabel<-as.vector(rowcol)
  res$ncorr<-0
  for (j in 1:min(k_clust,kclusttrue)) {
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



#' @title get_clusters
#'
#' @description Network-based clustering
#'
#' @param myData Data to be clustered, must be either binary (with levels "0"/"1") or categorical (with levels "0"/"1"/"2"/...)
#' @param k_clust Number of clusters
#' @param n_bg Number of covariates to be adjusted for; the position of the covariates must be in the last column of the myData matrix
#' @param EMseeds Seeds
#' @param edgepmat Matrix of penalized edges in the search space
#' @param blacklist Matrix of forbidden edges in the search space
#' @param bdepar Hyperparameters for structure learning (BDE score)
#' @param newallrelativeprobabs relative probability of cluster assignment of each sample
#' @param quick if TRUE, then the runtime is quick but accuracy is lower
# #' @param err error threshold defining when to stop
#'
#' @return a list containing the clusterMemberships and "assignprogress"
#' @export
#' @examples
#' \donttest{
#' # choose data
#' sampled_data <- sampleData(n_vars = 15, n_samples = c(300,300,300))$sampled_data
#' # learn clusters
#' cluster_results <- get_clusters(sampled_data)
#' # visualize the networks
#' library(ggplot2)
#' library(ggraph)
#' library(igraph)
#' library(ggpubr)
#' plot_clusters(cluster_results)
#' }
get_clusters <- function(myData, k_clust=3, n_bg=0, quick=TRUE, EMseeds=1, edgepmat=NULL, blacklist=NULL, bdepar=list(chi = 0.5, edgepf = 8), newallrelativeprobabs=NULL){

  # measure time
  start_time <- Sys.time()

  if (missing(myData)) stop("Need a categorical matrix as input to cluster.")
  if (!all(myData%%1==0)) stop("All categorical variables need to be specified as integers. Binary variables can be 0 or 1.")

  if (!all(myData < 2)){
    # cetegorical version
    score_type <- "bdecat"
  }else{
    # binary version
    score_type <- "bde"
  }

  #define when EM converges
  if(quick){
    err<-1e-6
    nIterations <- 10
    itLim <- 50
  }else{
    err<-1e-10
    nIterations <- 30
    itLim <- 100
    EMseeds <- EMseeds[1]+0:2
  }

  startseed <- EMseeds[1]

  # Prior Bernoulli clustering
  if(is.null(newallrelativeprobabs)){
    # nIterations <- 30
    chi <- 0.5 # pseudocounts for the Beta prior
    binClust <- get_clusters_bernoulli(binaryMatrix = myData, chi = chi, k_clust = k_clust, startseed = startseed, nIterations = nIterations, verbose=TRUE)
    newallrelativeprobabs <- binClust$relativeweights
    newclustermembership <- binClust$newclustermembership
    # newclustermembership <- reassignsamples(newallrelativeprobabs)
  }

  # # pre-clustering step
  # if (categorical){
  #   score_type <- "bdecat"
  #   # Categorical clustering
  #   if(is.null(newallrelativeprobabs)){
  #     nIterations <- 10
  #     chi <- 0.5 # pseudocounts for the Beta prior
  #     binClust <- BMMclusterEM(binaryMatrix = myData, chi = chi, k_clust = k_clust, startseed = startseed, nIterations = nIterations, verbose=TRUE)
  #     newallrelativeprobabs <- binClust$relativeweights
  #     newclustermembership <- binClust$newclustermembership
  #     newclustermembership <- reassignsamples(newallrelativeprobabs)
  #   }
  # }else{
  #   score_type <- "bde"
  #   # Binary clustering
  #   if(is.null(newallrelativeprobabs)){
  #     nIterations <- 10
  #     chi <- 0.5 # pseudocounts for the Beta prior
  #     binClust <- BBMMclusterEM(binaryMatrix = myData, chi = chi, k_clust = k_clust, startseed = startseed, nIterations = nIterations, verbose=TRUE)
  #     newallrelativeprobabs <- binClust$relativeweights
  #     newclustermembership <- binClust$newclustermembership
  #     newclustermembership <- reassignsamples(newallrelativeprobabs)
  #   }
  # }

  # # number of iterations (and respective seeds)
  # nSeeds <- 2
  # # number of covariate variables
  # n_bg <- 3
  ss <- dim(myData)[1]
  # #number of clusters
  # k_clust<-3
  # # iteration limit
  # itLim <- 20

  # total number of variables
  nn <- dim(myData)[2]
  #number of variables (without covariates)
  n<-nn-n_bg

  # #prior pseudo counts
  # chixi<-0.5
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
  #to store memberships
  newclustermembership <- list()
  #to store scores against clusters
  scoresagainstclusters<-matrix(ncol=k_clust,nrow=ss)
  # initial cluster assignment
  initial_newallrelativeprobabs <- newallrelativeprobabs
  # store clustercenters
  all_clustercenters <- list()

  for (s in 1:EMn) {
    diffy<-1
    cnt<-1
    newallrelativeprobabs <- initial_newallrelativeprobabs
    assignprogress[[s]]<-list()
    assignprogress_local<-list()
    # assignprogress_local$corrvec<-numeric()
    assignprogress_local$likel<-numeric()
    # set.seed(EMseeds[s])
    print(paste("EM seed",EMseeds[s]))
    # centers<-list()
    for (i in 1:k_clust) {
      # centers[[i]]<-as.data.frame(matrix(ncol=3))
      clustercenters[[i]]<-matrix(rep(0,(n+n_bg)^2),nrow=n+n_bg)
    }

    # if(!BBMMClust){
    #   #generate random assignment of belonging to each cluster for each sample
    #   newallrelativeprobabs<-generatetriple(ss)
    #   newclustermembership<-reassignsamples(newallrelativeprobabs)
    # }


    # #learn how many samples were asssigned correctly
    # res<-checkmembership(k_clust,kclusttrue,truelabels,newclustermembership)
    # print(paste("number of correctly assigned samples by random assignment:", res$ncorr,
    #             "of total", ss, "samples"))

    # EM cycle
    while (diffy>err&cnt<itLim) {
      allrelativeprobabs<-newallrelativeprobabs
      allprobprev<-newallrelativeprobabs
      coltots<-colSums(allrelativeprobabs) + bdepar$chi # add prior to clustersizes
      tauvec<-coltots/sum(coltots)

      #outer EM cycle, learning cluster centers
      #define cluster centers and assign probabilities

      clustercenters <- parallel::mclapply(1:k_clust, function(k) {

        #define score parameters
        if (n_bg>0){
          # apply different scores for categorical / binary data
          scorepar <- BiDAG::scoreparameters(score_type, as.data.frame(myData), edgepmat = edgepmat,
                                           weightvector=allrelativeprobabs[,k],
                                           bdepar=bdepar, bgnodes=(n+1):(n+n_bg))
        }else{
          scorepar <- BiDAG::scoreparameters(score_type, as.data.frame(myData), edgepmat = edgepmat,
                                           weightvector=allrelativeprobabs[,k],
                                           bdepar=bdepar)
        }

        # check mark

        #find MAP DAG using iterative order MCMC
        maxfit<-BiDAG::iterativeMCMC(scorepar,addspace=clustercenters[[k]],verbose=FALSE,blacklist=blacklist)
        #store it
        clustercenters[[k]]<-maxfit$DAG

        return(clustercenters[[k]])

      # }, mc.cores = k_clust)
      })

      #internal EM cycle, estimating parameters
      for (i in 1:nit.internal) {
        allrelativeprobabs<-newallrelativeprobabs
        coltots<-colSums(allrelativeprobabs) + bdepar$chi # add prior to clustersizes
        tauvec<-coltots/sum(coltots)

        parRes <- parallel::mclapply(1:k_clust, function(k) {
          if (n_bg>0){
            scorepar <- BiDAG::scoreparameters(score_type,as.data.frame(myData), edgepmat = edgepmat,
                                             weightvector=allrelativeprobabs[,k],
                                             bdepar=bdepar, bgnodes=(n+1):(n+n_bg))
          }else{
            scorepar <- BiDAG::scoreparameters(score_type,as.data.frame(myData), edgepmat = edgepmat,
                                             weightvector=allrelativeprobabs[,k],
                                             bdepar=bdepar)
          }

          scorepar$n <- n # to avoid to scoring over background nodes
          # scoresagainstclusters[,k] <- BiDAG::scoreagainstDAGscoreagainstDAG(scorepar,clustercenters[[k]])

          if (score_type=="bdecat"){
            scoresagainstclusters[,k] <- BiDAG::scoreagainstDAG(scorepar,clustercenters[[k]], bdecatCvec = apply(myData, 2, function(x) length(unique(x))))
          }else{
            scoresagainstclusters[,k] <- BiDAG::scoreagainstDAG(scorepar,clustercenters[[k]])
          }

          scorepar$n <- n+n_bg # recet after scoring

          return(scoresagainstclusters[,k])

        # }, mc.cores = k_clust)
        })

        for (kk in 1:k_clust){
          scoresagainstclusters[,kk] <- parRes[[kk]]
        }

        newallrelativeprobabsnotau <- allrelativeprobs(scoresagainstclusters)
        newallrelativeprobabs <- relativeprobswithtau(newallrelativeprobabsnotau,tauvec)
      }

      diffy<-sum((allprobprev-newallrelativeprobabs)^2)
      newclustermembership[[s]]<-reassignsamples(newallrelativeprobabs)
      # res<-checkmembership(k_clust,kclusttrue=3,truelabels,newclustermembership)
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
    all_clustercenters[[s]] <- clustercenters
  }

  # assignprogress[[s]]$likel[length(cluster_res_t$assignprogress$likel)]

  if(EMn>1){
    bestSeed <- getBestSeed(assignprogress)
    newallrelativeprobabs <- probs[[bestSeed]]
    assignprogress <- assignprogress[[bestSeed]][[1]]
    clustercenters <- all_clustercenters[[bestSeed]]
    newclustermembership <- newclustermembership[[bestSeed]]
    names(newclustermembership) <- rownames(myData)
    print(paste0("Best seed: ", bestSeed))
  }else{
    newclustermembership <- newclustermembership[[1]]
    assignprogress <- assignprogress[[1]]
  }

  #   # get best performing seed
  #   assignprogressList <- lapply(clusterResAll, function(x) x[[2]][[1]])
  #   bestSeed <- getBestSeed(assignprogressList)
  #   bestRes <- clusterResAll[[bestSeed]]

  # measure time
  end_time <- Sys.time()

  print(paste0("Computation time: ",round(as.numeric(difftime(end_time, start_time, units='mins')), digits = 2), " mins"))

  return(list("clustermembership"=newclustermembership,"assignprogress"=assignprogress, "DAGs"=clustercenters, "probs"=newallrelativeprobabs, "data"=myData, "n_bg"=n_bg))
}



# #' @title get_clusters_binary
# #'
# #' @description Network-based clustering
# #'
# #' @param myData Data to be clustered
# #' @param k_clust Number of clusters
# #' @param n_bg Number of covariates
# #' @param itLim Maximum number of iterations
# #' @param EMseeds Seeds
# #' @param edgepmat Matrix of penalized edges in the search space
# #' @param bdepar Hyperparameters for structure learning (BDE score)
# #' @param newallrelativeprobabs relative probability of cluster assignment of each sample
# #'
# #' @return a list containing the clusterMemberships and "assignprogress"
# #' @export
# get_clusters_binary <- function(myData,k_clust=3,n_bg=0,itLim=50, EMseeds=1, edgepmat=NULL, bdepar=list(chi = 0.5, edgepf = 16), newallrelativeprobabs=NULL){
#
#   # measure time
#   start_time <- Sys.time()
#
#   # Binary clustering
#   startseed <- EMseeds[1]
#   nIterations <- 10
#   chi <- 0.5 # pseudocounts for the Beta prior
#
#   if(is.null(newallrelativeprobabs)){
#     binClust <- BBMMclusterEM(binaryMatrix = myData, chi = chi, k_clust = k_clust, startseed = startseed, nIterations = nIterations, verbose=TRUE)
#     newallrelativeprobabs <- binClust$relativeweights
#     newclustermembership <- binClust$newclustermembership
#     newclustermembership <- reassignsamples(newallrelativeprobabs)
#   }
#
#   # # number of iterations (and respective seeds)
#   # nSeeds <- 2
#   # # number of covariate variables
#   # n_bg <- 3
#   ss <- dim(myData)[1]
#   # #number of clusters
#   # k_clust<-3
#   # # iteration limit
#   # itLim <- 20
#
#   # total number of variables
#   nn <- dim(myData)[2]
#   #number of variables (without covariates)
#   n<-nn-n_bg
#
#   # #prior pseudo counts
#   # chixi<-0.5
#   #define when EM converges
#   err<-1e-6
#   # #edge penalization factor
#   # edgepfy<-16
#   #define different seeds to run EM
#   # EMseeds<-c(101,102,103,104,105)
#   # EMseeds<-c(100)+c(1:nSeeds)
#   #number of EM attempts to ensure highest likelihood
#   EMn<-length(EMseeds)
#   #number of iterations in the internal cycle
#   nit.internal<-10
#
#
#   #to store accuracy of cluster centers
#   # centerprogress<-list()
#   #accuracy of assignments for different EM runs
#   assignprogress<-list()
#   #to store scores against clusters for each EM run
#   scoresprogress<-list()
#   #to store cluster probabilities for all sample for each EM run
#   probs<-list()
#   #to store relabelling
#   relabs<-list()
#   #to store cluster centers
#   clustercenters<-list()
#   #to store scores against clusters
#   scoresagainstclusters<-matrix(ncol=k_clust,nrow=ss)
#
#   for (s in 1:EMn) {
#     diffy<-1
#     cnt<-1
#     assignprogress[[s]]<-list()
#     assignprogress_local<-list()
#     # assignprogress_local$corrvec<-numeric()
#     assignprogress_local$likel<-numeric()
#     # set.seed(EMseeds[s])
#     print(paste("EM seed",EMseeds[s]))
#     # centers<-list()
#     for (i in 1:k_clust) {
#       # centers[[i]]<-as.data.frame(matrix(ncol=3))
#       clustercenters[[i]]<-matrix(rep(0,(n+n_bg)^2),nrow=n+n_bg)
#     }
#
#
#     # if(!BBMMClust){
#     #   #generate random assignment of belonging to each cluster for each sample
#     #   newallrelativeprobabs<-generatetriple(ss)
#     #   newclustermembership<-reassignsamples(newallrelativeprobabs)
#     # }
#
#
#     # #learn how many samples were asssigned correctly
#     # res<-checkmembership(k_clust,kclusttrue,truelabels,newclustermembership)
#     # print(paste("number of correctly assigned samples by random assignment:", res$ncorr,
#     #             "of total", ss, "samples"))
#
#     # EM cycle
#     while (diffy>err&cnt<itLim) {
#       allrelativeprobabs<-newallrelativeprobabs
#       allprobprev<-newallrelativeprobabs
#       coltots<-colSums(allrelativeprobabs) + bdepar$chi # add prior to clustersizes
#       tauvec<-coltots/sum(coltots)
#
#       #outer EM cycle, learning cluster centers
#       #define cluster centers and assign probabilities
#
#       clustercenters <- parallel::mclapply(1:k_clust, function(k) {
#
#         #define score parameters
#         if (n_bg>0){
#           # apply different scores for categorical / binary data
#             scorepar <- BiDAG::scoreparameters("bde",myData, edgepmat = edgepmat,
#                                                weightvector=allrelativeprobabs[,k],
#                                                bdepar=bdepar, bgnodes=(n+1):(n+n_bg))
#         }else{
#             scorepar <- BiDAG::scoreparameters("bde",myData, edgepmat = edgepmat,
#                                                weightvector=allrelativeprobabs[,k],
#                                                bdepar=bdepar)
#         }
#
#         #find MAP DAG using iterative order MCMC
#         maxfit<-BiDAG::iterativeMCMC(scorepar,addspace=clustercenters[[k]],verbose=FALSE,blacklist=blacklist)
#         #store it
#         clustercenters[[k]]<-maxfit$DAG
#
#         return(clustercenters[[k]])
#
#       }, mc.cores = k_clust)
#
#       #internal EM cycle, estimating parameters
#       for (i in 1:nit.internal) {
#         allrelativeprobabs<-newallrelativeprobabs
#         coltots<-colSums(allrelativeprobabs) + bdepar$chi # add prior to clustersizes
#         tauvec<-coltots/sum(coltots)
#
#         parRes <- parallel::mclapply(1:k_clust, function(k) {
#           if (n_bg>0){
#             scorepar <- BiDAG::scoreparameters("bde",as.data.frame(myData), edgepmat = edgepmat,
#                                                weightvector=allrelativeprobabs[,k],
#                                                bdepar=bdepar, bgnodes=(n+1):(n+n_bg))
#           }else{
#             scorepar <- BiDAG::scoreparameters("bde",as.data.frame(myData), edgepmat = edgepmat,
#                                                weightvector=allrelativeprobabs[,k],
#                                                bdepar=bdepar)
#           }
#           scorepar$n <- n # to avoid to scoring over background nodes
#           scoresagainstclusters[,k] <- BiDAG::scoreagainstDAG(scorepar,clustercenters[[k]])
#           scorepar$n <- n+n_bg # recet after scoring
#
#           return(scoresagainstclusters[,k])
#
#         }, mc.cores = k_clust)
#
#         for (kk in 1:k_clust){
#           scoresagainstclusters[,kk] <- parRes[[kk]]
#         }
#
#         newallrelativeprobabsnotau<-allrelativeprobs(scoresagainstclusters)
#         newallrelativeprobabs<-relativeprobswithtau(newallrelativeprobabsnotau,tauvec)
#       }
#
#       diffy<-sum((allprobprev-newallrelativeprobabs)^2)
#       newclustermembership<-reassignsamples(newallrelativeprobabs)
#       # res<-checkmembership(k_clust,kclusttrue=3,truelabels,newclustermembership)
#       # assignprogress_local$corrvec[cnt]<-res$ncorr
#       assignprogress_local$likel[cnt]<-calcloglike(scoresagainstclusters,tauvec)
#       # for (k in 1:3) {
#       #   centers[[k]][cnt,]<-compareDAGs(m2graph(clustercenters[[k]]),
#       #                                   BNsBG[[res$relabel[k]]]$DAG)[c("SHD","TP","FP")]
#       # }
#       cnt<-cnt+1
#     }
#     print(paste("EM converged after",cnt-1,"iterations"))
#     # print(paste("number of correctly assigned samples:",
#     # res$ncorr, "of", ss))
#     #store number of correct assignments and likelihood
#     assignprogress[[s]]<-assignprogress_local
#     #store center progress
#     # centerprogress[[s]]<-centers
#     scoresprogress[[s]]<-scoresagainstclusters
#     probs[[s]]<-newallrelativeprobabs
#     # relabs[[s]]<-res$relabel
#   }
#
#   # measure time
#   end_time <- Sys.time()
#
#   print(paste0("Computation time: ",as.numeric(end_time-start_time)))
#
#   names(newclustermembership) <- rownames(myData)
#
#   return(list("clustermembership"=newclustermembership,"assignprogress"=assignprogress, "DAGs"=clustercenters, "newallrelativeprobabs"=newallrelativeprobabs))
# }


# #' @title clustNetParallel
# #'
# #' @description Network-based clustering of multiple seeds using parallel computing
# #'
# #' @param myData Data to be clustered
# #' @param k_clust Number of clusters
# #' @param n_bg Number of covariates
# #' @param itLim Maximum number of iterations
# #' @param EMseeds seeds
# #' @param edgepmat a matrix of penalized edges in the search space
# #' @param bdepar Hyperparameters for structure learning (BDE score)
# #'
# #' @return a list containing the clusterMemberships, DAGs, best seed and "assignprogress"
# #' @export
# #'
# clustNetParallel <- function(myData,k_clust=3,n_bg=0,itLim=20, EMseeds=1:5, edgepmat=NULL, bdepar=list(chi = 0.5, edgepf = 16)){
#
#   # parallel computing of clustering
#   nSeeds <- length(EMseeds)
#   clusterResAll <- parallel::mclapply(1:nSeeds, function(i) {
#     print(paste("Clustering iteration", i, "of", nSeeds))
#     clusterRes <- get_clusters(myData=myData,k_clust=k_clust,n_bg=n_bg,itLim=itLim, EMseeds=EMseeds[i], edgepmat=edgepmat, bdepar=bdepar)
#     return(clusterRes)
#   }, mc.cores = nSeeds)
#
#   # get best performing seed
#   assignprogressList <- lapply(clusterResAll, function(x) x[[2]][[1]])
#   bestSeed <- getBestSeed(assignprogressList)
#   bestRes <- clusterResAll[[bestSeed]]
#
#   # return results of best seed
#   return(list("clustermembership"=bestRes$clustermembership,"assignprogress"=bestRes$assignprogress,"DAGs"=bestRes$DAGs,"bestSeed"=bestSeed))
# }

# #' @title plot_clusters
# #'
# #' @description Plot DAGs of each cluster
# #'
# #' @param learned_clusters Output from get_clusters()
# #'
# #' @return A plot of DAGs
# #' @export
# #'
# plot_clusters <- function(learned_clusters){
#   # plot DAGs of each cluster
#   n_dags <- length(learned_clusters$DAGs)
#   par(mfrow = c(2, 1+as.integer(n_dags/2)))
#   for (ii in 1:n_dags){
#     SGS:::plot_dag(as.matrix(learned_clusters$DAGs[ii][[1]]))
#   }
#   par(mfrow = c(1,1))
# }


# #' @title clustNeter
# #'
# #' @description Network-based clustering
# #'
# #' @param myData Data to be clustered
# #' @param k_clust Number of clusters
# #' @param n_bg Number of covariates
# #' @param itLim Maximum number of iterations
# #' @param EMseeds Seeds
# #' @param BBMMClust Binary clustering before network-based clustering (TRUE by default)
# #' @param edgepmat Matrix of penalized edges in the search space
# #' @param bdepar Hyperparameters for structure learning (BDE score)
# #'
# #' @return a list containing the clusterMemberships and "assignprogress"
# #' @export
# clustNeter <- function(myData,k_clust=3,n_bg=0,itLim=20, EMseeds=1, BBMMClust=TRUE, edgepmat=NULL, bdepar=list(chi = 0.5, edgepf = 16)){
#
#   # Binary clustering
#   startseed <- EMseeds[1]
#   nIterations <- 10
#   chi <- 0.5 # pseudocounts for the Beta prior
#
#   if(BBMMClust){
#     binClust <- BBMMclusterEM(binaryMatrix = myData, chi = chi, k_clust = k_clust, startseed = startseed, nIterations = nIterations, verbose=TRUE)
#     newallrelativeprobabs <- binClust$relativeweights
#     newclustermembership <- binClust$newclustermembership
#     newclustermembership <- reassignsamples(newallrelativeprobabs)
#   }
#
#   # # number of iterations (and respective seeds)
#   # nSeeds <- 2
#   # # number of covariate variables
#   # n_bg <- 3
#   ss <- dim(myData)[1]
#   # #number of clusters
#   # k_clust<-3
#   # # iteration limit
#   # itLim <- 20
#
#   # total number of variables
#   nn <- dim(myData)[2]
#   #number of variables (without covariates)
#   n<-nn-n_bg
#
#   # #prior pseudo counts
#   # chixi<-0.5
#   #define when EM converges
#   err<-1e-6
#   # #edge penalization factor
#   # edgepfy<-16
#   #define different seeds to run EM
#   # EMseeds<-c(101,102,103,104,105)
#   # EMseeds<-c(100)+c(1:nSeeds)
#   #number of EM attempts to ensure highest likelihood
#   EMn<-length(EMseeds)
#   #number of iterations in the internal cycle
#   nit.internal<-10
#
#
#   #to store accuracy of cluster centers
#   # centerprogress<-list()
#   #accuracy of assignments for different EM runs
#   assignprogress<-list()
#   #to store scores against clusters for each EM run
#   scoresprogress<-list()
#   #to store cluster probabilities for all sample for each EM run
#   probs<-list()
#   #to store relabelling
#   relabs<-list()
#   #to store cluster centers
#   clustercenters<-list()
#   #to store scores against clusters
#   scoresagainstclusters<-matrix(ncol=k_clust,nrow=ss)
#
#   for (s in 1:EMn) {
#     diffy<-1
#     cnt<-1
#     assignprogress[[s]]<-list()
#     assignprogress_local<-list()
#     # assignprogress_local$corrvec<-numeric()
#     assignprogress_local$likel<-numeric()
#     set.seed(EMseeds[s])
#     print(paste("EM seed",EMseeds[s]))
#     # centers<-list()
#     for (i in 1:k_clust) {
#       # centers[[i]]<-as.data.frame(matrix(ncol=3))
#       clustercenters[[i]]<-matrix(rep(0,(n+n_bg)^2),nrow=n+n_bg)
#     }
#
#
#     if(!BBMMClust){
#       #generate random assignment of belonging to each cluster for each sample
#       newallrelativeprobabs<-generatetriple(ss)
#       newclustermembership<-reassignsamples(newallrelativeprobabs)
#     }
#
#
#     # #learn how many samples were asssigned correctly
#     # res<-checkmembership(k_clust,kclusttrue,truelabels,newclustermembership)
#     # print(paste("number of correctly assigned samples by random assignment:", res$ncorr,
#     #             "of total", ss, "samples"))
#
#     # EM cycle
#     while (diffy>err&cnt<itLim) {
#       allrelativeprobabs<-newallrelativeprobabs
#       allprobprev<-newallrelativeprobabs
#       coltots<-colSums(allrelativeprobabs) + bdepar$chi # add prior to clustersizes
#       tauvec<-coltots/sum(coltots)
#
#       #outer EM cycle, learning cluster centers
#       #define cluster centers and assign probabilities
#       for (k in 1:k_clust) {
#
#         # #learn background nodes
#         # for (i_n_bg in 1:n_bg){
#         #   for (i_n in 1:n){
#         #     #chi squared test
#         #     testRes <- wtd.chi.sq(myData[,i_n],myData[,n+i_n_bg],weight=allrelativeprobabs[,k])
#         #     siglevel <- 0.05/(k_clust+n)
#         #     #add to initial DAG structure
#         #     if (testRes[3]<siglevel){
#         #       clustercenters[[k]][n+i_n_bg,i_n] <- 1
#         #     }else{
#         #       clustercenters[[k]][n+i_n_bg,i_n] <- 0
#         #     }
#         #   }
#         # }
#
#         #define score parameters
#         if (n_bg>0){
#           scorepar<-BiDAG::scoreparameters("bde",myData, edgepmat = edgepmat,
#                                     weightvector=allrelativeprobabs[,k],
#                                     bdepar=bdepar, bgnodes=(n+1):(n+n_bg))
#         }else{
#           scorepar<-BiDAG::scoreparameters("bde",myData, edgepmat = edgepmat,
#                                     weightvector=allrelativeprobabs[,k],
#                                     bdepar=bdepar)
#         }
#
#         #find MAP DAG using iterative order MCMC
#         maxfit<-BiDAG::iterativeMCMC(scorepar,addspace=clustercenters[[k]],verbose=FALSE,blacklist=blacklist)
#         #store it
#         clustercenters[[k]]<-maxfit$DAG#
#       }
#
#       #internal EM cycle, estimating parameters
#       for (i in 1:nit.internal) {
#         allrelativeprobabs<-newallrelativeprobabs
#         coltots<-colSums(allrelativeprobabs) + bdepar$chi # add prior to clustersizes
#         tauvec<-coltots/sum(coltots)
#         for (k in 1:k_clust) {
#           if (n_bg>0){
#             scorepar<-BiDAG::scoreparameters("bde",as.data.frame(myData), edgepmat = edgepmat,
#                                              weightvector=allrelativeprobabs[,k],
#                                              bdepar=bdepar, bgnodes=(n+1):(n+n_bg))
#           }else{
#             scorepar<-BiDAG::scoreparameters("bde",as.data.frame(myData), edgepmat = edgepmat,
#                                              weightvector=allrelativeprobabs[,k],
#                                              bdepar=bdepar)
#           }
#           scorepar$n <- n # to avoid to scoring over background nodes
#           scoresagainstclusters[,k]<-BiDAG::scoreagainstDAG(scorepar,clustercenters[[k]])
#           scorepar$n <- n+n_bg # recet after scoring
#         }
#         newallrelativeprobabsnotau<-allrelativeprobs(scoresagainstclusters)
#         newallrelativeprobabs<-relativeprobswithtau(newallrelativeprobabsnotau,tauvec)
#       }
#
#       diffy<-sum((allprobprev-newallrelativeprobabs)^2)
#       newclustermembership<-reassignsamples(newallrelativeprobabs)
#       # res<-checkmembership(k_clust,kclusttrue=3,truelabels,newclustermembership)
#       # assignprogress_local$corrvec[cnt]<-res$ncorr
#       assignprogress_local$likel[cnt]<-calcloglike(scoresagainstclusters,tauvec)
#       # for (k in 1:3) {
#       #   centers[[k]][cnt,]<-compareDAGs(m2graph(clustercenters[[k]]),
#       #                                   BNsBG[[res$relabel[k]]]$DAG)[c("SHD","TP","FP")]
#       # }
#       cnt<-cnt+1
#     }
#     print(paste("EM converged after",cnt-1,"iterations"))
#     # print(paste("number of correctly assigned samples:",
#     # res$ncorr, "of", ss))
#     #store number of correct assignments and likelihood
#     assignprogress[[s]]<-assignprogress_local
#     #store center progress
#     # centerprogress[[s]]<-centers
#     scoresprogress[[s]]<-scoresagainstclusters
#     probs[[s]]<-newallrelativeprobabs
#     # relabs[[s]]<-res$relabel
#   }
#
#   return(list("clustermembership"=newclustermembership,"assignprogress"=assignprogress, "DAGs"=clustercenters))
# }


# #' @title clustNeterParallel
# #'
# #' @description Network-based clustering of multiple seeds using parallel computing
# #'
# #' @param myData Data to be clustered
# #' @param k_clust Number of clusters
# #' @param n_bg Number of covariates
# #' @param itLim Maximum number of iterations
# #' @param EMseeds seeds
# #' @param BBMMClust binary clustering before network-based clustering (TRUE by default)
# #' @param edgepmat a matrix of penalized edges in the search space
# #' @param bdepar Hyperparameters for structure learning (BDE score)
# #'
# #' @return a list containing the clusterMemberships, DAGs, best seed and "assignprogress"
# #' @export
# #'
# clustNeterParallel <- function(myData,k_clust=3,n_bg=0,itLim=20, EMseeds=1:5, BBMMClust=TRUE, edgepmat=NULL, bdepar=list(chi = 0.5, edgepf = 16)){
#
#   # parallel computing of clustering
#   nSeeds <- length(EMseeds)
#   clusterResAll <- parallel::mclapply(1:nSeeds, function(i) {
#     print(paste("Clustering iteration", i, "of", nSeeds))
#     clusterRes <- clustNeter(myData=myData,k_clust=k_clust,n_bg=n_bg,itLim=itLim, EMseeds=EMseeds[i], BBMMClust=BBMMClust, edgepmat=edgepmat, bdepar=bdepar)
#     return(clusterRes)
#   }, mc.cores = nSeeds)
#
#   # get best performing seed
#   assignprogressList <- lapply(clusterResAll, function(x) x[[2]][[1]])
#   bestSeed <- getBestSeed(assignprogressList)
#   bestRes <- clusterResAll[[bestSeed]]
#
#   # return results of best seed
#   return(list("clustermembership"=bestRes$clustermembership,"assignprogress"=bestRes$assignprogress,"DAGs"=bestRes$DAGs,"bestSeed"=bestSeed))
# }

