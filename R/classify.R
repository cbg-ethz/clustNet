#' @title get_classification
#'
#' @description Classification based on clustering
#'
#' @param cluster_results Output from get_clusters()
#' @param data_classify Data that should be classified; colnames need to match the ones of cluster_results$data; missing cols are allowed
#'
#' @return a list containing the classification as "clustermembership" and the probabilities of belonging to the clusters as "allrelativeprobabs"
#' @export
#' @examples
#' \donttest{
#' # choose data
#' sampled_data <- sampleData(n_vars = 15, n_samples = c(300,300,300))$sampled_data
#' # learn clusters
#' cluster_results <- get_clusters(sampled_data)
#' # visualize the networks
#' classification_results <- get_classification(cluster_results, sampled_data)
#' }
get_classification <- function(cluster_results, data_classify){

  myData <- cluster_results$data
  k_clust <- length(cluster_results$DAGs)

  if(is.vector(data_classify)){
    data_classify <- t(as.data.frame(data_classify)) # when this is a single col entry
  }else{
    data_classify <- as.data.frame(data_classify)
  }

  # input is clustercenters
  clustercenters <- cluster_results$DAGs
  newallrelativeprobabs <- cluster_results$probs


  ## detect and adjust for missing data

  # Create two example matrices
  matrix_one <- myData

  # Create two example matrices
  matrix_two <- data_classify

  # Find shared variables
  shared_vars <- intersect(colnames(matrix_one), colnames(matrix_two))
  # cat("Shared variables:", shared_vars, "\n")

  # Find missing variables in matrix_one
  missing_vars_one <- setdiff(colnames(matrix_one), colnames(matrix_two))
  if(length(missing_vars_one)> 0){
    cat("Missing variables in cluster data:", missing_vars_one, "\n")
  }

  # Find missing variables in matrix_two
  missing_vars_two <- setdiff(colnames(matrix_two), colnames(matrix_one))
  if(length(missing_vars_two)> 0){
    cat("Missing variables in classification data:", missing_vars_two, "\n")
  }

  # Get column indices of shared variables in matrix_one
  if (length(shared_vars) > 0) {
    shared_columns_indices_one <- match(shared_vars, colnames(matrix_one))
    shared_columns_indices_two <- match(shared_vars, colnames(matrix_two))
    # adapt data to shared vars
    myData <- matrix_one[,shared_columns_indices_one]
    data_classify <- matrix_two[,shared_columns_indices_two]
    # adapt DAGs to shared vars
    for (ii in 1:k_clust){
      clustercenters[[ii]] <- clustercenters[[ii]][shared_columns_indices_one,shared_columns_indices_one]
    }
  } else {
    cat("No shared variables\n")
  }

  ## classification

  data_classify <- unname(data_classify)

  # number of background variables
  n_bg <- 2
  # total number of variables
  ss <- dim(myData)[1]
  nn <- dim(myData)[2]
  ss2 <- NROW(data_classify)
  #number of variables (without covariates)
  n<-nn-n_bg

  scoresagainstclusters<-matrix(ncol=k_clust,nrow=ss2)

  bdepar <- list(chi = 0.5, edgepf = 8)


  if (!all(myData < 2)){
    # cetegorical version
    score_type <- "bdecat"
  }else{
    # binary version
    score_type <- "bde"
  }

  allrelativeprobabs<-newallrelativeprobabs
  coltots<-colSums(allrelativeprobabs) + bdepar$chi # add prior to clustersizes
  tauvec<-coltots/sum(coltots)

  parRes <- parallel::mclapply(1:k_clust, function(k) {
    if (n_bg>0){
      scorepar <- BiDAG::scoreparameters(score_type,as.data.frame(myData),
                                         weightvector=allrelativeprobabs[,k],
                                         bdepar=bdepar, bgnodes=(n+1):(n+n_bg))
    }else{
      scorepar <- BiDAG::scoreparameters(score_type,as.data.frame(myData),
                                         weightvector=allrelativeprobabs[,k],
                                         bdepar=bdepar)
    }

    scorepar$n <- n # to avoid to scoring over background nodes
    # scoresagainstclusters[,k] <- BiDAG::scoreagainstDAGscoreagainstDAG(scorepar,clustercenters[[k]])

    # if (score_type=="bdecat"){
    #   scoresagainstclusters[,k] <- BiDAG::scoreagainstDAG(scorepar,clustercenters[[k]], bdecatCvec = apply(myData, 2, function(x) length(unique(x))))
    # }else{
    #   scoresagainstclusters[,k] <- BiDAG::scoreagainstDAG(scorepar,clustercenters[[k]])
    # }

    if (score_type=="bdecat"){
      scoresagainstclusters[,k] <- BiDAG::scoreagainstDAG(scorepar,clustercenters[[k]], bdecatCvec = apply(myData, 2, function(x) length(unique(x))), datatoscore = data_classify)
    }else{
      scoresagainstclusters[,k] <- BiDAG::scoreagainstDAG(scorepar,clustercenters[[k]], datatoscore = data_classify)
    }


    scorepar$n <- n+n_bg # recet after scoring

    return(scoresagainstclusters[,k])

    # }, mc.cores = k_clust)
  })

  for (kk in 1:k_clust){
    scoresagainstclusters[,kk] <- parRes[[kk]]
  }

  # assign cluster
  newallrelativeprobabsnotau <- allrelativeprobs(scoresagainstclusters)
  newallrelativeprobabs <- relativeprobswithtau(newallrelativeprobabsnotau,tauvec)

  newclustermembership<-reassignsamples(newallrelativeprobabs)

  return(list("clustermembership"=newclustermembership,"allrelativeprobabs"=newallrelativeprobabs, "shared_vars"=shared_vars, "classified_data"=myData))
}
