function (scorepar, incidence, datatoscore = NULL, marginalise = FALSE, 
          onlymain = FALSE, bdecatCvec = NULL) 
{
  if (onlymain & !marginalise) {
    n <- scorepar$n
    mainnodes <- scorepar$mainnodes
  }
  else {
    n <- scorepar$n
    mainnodes <- c(1:n)
  }
  if (is.null(datatoscore)) {
    datatoscore <- scorepar$data
  }
  if (is.vector(datatoscore)) {
    datatoscore <- matrix(datatoscore, nrow = 1)
  }
  if (scorepar$type == "bge") {
    if (marginalise == FALSE) {
      datatoscore <- t(t(datatoscore) - scorepar$muN)
    }
    else {
      datatoscore <- t(t(datatoscore) - scorepar$means)
    }
  }
  if (scorepar$type == "bge" && marginalise != FALSE) {
    return(scoreagainstDAGmargBGe(n, scorepar, incidence, 
                                  datatoscore))
  }
  else if (scorepar$type == "mixed") {
    binscore <- scoreagainstDAG(scorepar$binpar, incidence[1:scorepar$nbin, 
                                                           1:scorepar$nbin])
    gausscore <- scoreagainstDAG(scorepar$gausspar, incidence)
    return(binscore + gausscore)
  }
  else if (scorepar$type == "bde") {
    samplescores <- matrix(0, nrow = nrow(datatoscore), ncol = n)
    for (j in mainnodes) {
      parentnodes <- which(incidence[, j] == 1)
      samplescores[, j] <- scoreagainstDAGcore(j, parentnodes, 
                                               n, scorepar, datatoscore)
    }
  }
  else {
    if (!is.null(bdecatCvec)) {
      scorepar$Cvec <- bdecatCvec
    }
    samplescores <- matrix(0, nrow = nrow(datatoscore), ncol = n)
    for (j in 1:n) {
      parentnodes <- which(incidence[, j] == 1)
      samplescores[, j] <- scoreagainstDAGcore(j, parentnodes, 
                                               n, scorepar, datatoscore)
    }
  }
  return(rowSums(samplescores))
}

function(scorepar, incidence, datatoscore=NULL, marginalise=FALSE, bdecatCvec=NULL){
  
  n<-scorepar$n
  
  if (is.null(datatoscore)) {
    datatoscore<-scorepar$data
  }
  
  if(is.vector(datatoscore)){ # if input is a vector
    datatoscore <- matrix(datatoscore, nrow=1) # cast it as a matrix
  }
  
  if (scorepar$type=="bge") {
    if(marginalise==FALSE){
      datatoscore <- t(t(datatoscore) - scorepar$muN) # recentre around posterior mean
    } else {
      datatoscore <- t(t(datatoscore) - scorepar$means) # recentre about data mean
    }
  }
  
  if (scorepar$type=="bge" && marginalise!=FALSE){
    return(BiDAG:::scoreagainstDAGmargBGe(n, scorepar, incidence, datatoscore))
  } else if (scorepar$type=="mixed") {
    binscore<-scoreagainstDAG(scorepar$binpar, incidence[1:scorepar$nbin,1:scorepar$nbin])
    gausscore<-scoreagainstDAG(scorepar$gausspar, incidence)
    return(binscore+gausscore)
    
  } else if (scorepar$type=="bde"){
    samplescores <- matrix(0,nrow=nrow(datatoscore),ncol=n)
    for (j in 1:n)  {
      parentnodes <- which(incidence[,j]==1)
      samplescores[,j]<-scoreagainstDAGcore(j,parentnodes,n,scorepar,datatoscore)
    }
    
    return(rowSums(samplescores))
  } else {
    if (!is.null(bdecatCvec)) {
      scorepar$Cvec <- bdecatCvec
    }
    samplescores <- matrix(0,nrow=nrow(datatoscore),ncol=n)
    for (j in 1:n)  {
      parentnodes <- which(incidence[,j]==1)
      samplescores[,j]<-scoreagainstDAGcore(j,parentnodes,n,scorepar,datatoscore)
    }
    
    return(rowSums(samplescores))
  }
}