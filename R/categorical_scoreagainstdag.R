#'Calculating the score of a sample against a DAG
#'
#'This function calculates the score of a given sample against a DAG represented by its incidence matrix.
#'@param scorepar an object of class \code{scoreparameters}; see constructor function \code{\link{scoreparameters}}
#'@param incidence a square matrix of dimensions equal to the number of variables with entries in \code{\{0,1\}}, representing the adjacency matrix of the DAG against which the score is calculated
#'@param datatoscore (optional) a matrix (vector) containing binary (for BDe score) or continuous (for the BGe score) observations (or just one observation)  to be scored; the number of columns should be equal to the number of variables in the Bayesian network, the number of rows should be equal to the number of observations; by default all data from \code{scorepar} parameter is used
#'@param marginalise (optional for continuous data), whether to use the posterior mean for scoring (default) or to marginalise over the posterior distribution (more computationally costly)
#'@return the log of the BDe/BGe score of given observations against a DAG
#'@examples
#'  Asiascore<-scoreparameters("bde", Asia[1:100,]) #we wish to score only first 100 observations
#'  scoreagainstDAG(Asiascore, Asiamat)
#'
#'@export

scoreagainstDAG <- function(scorepar, incidence, datatoscore=NULL, marginalise=FALSE, bdecatCvec=NULL){

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
    return(scoreagainstDAGmargBGe(n, scorepar, incidence, datatoscore))
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

# this function scores a nodes against its parents based on the BGe or BDe (binary) score


scoreagainstDAGcore<-function(j,parentnodes,n,param,datatoscore) {
  samplenodescores<-rep(0,nrow(datatoscore)) # store
  lp<-length(parentnodes) # number of parents

  switch(param$type,
         "bge" = {
           Sigma <- param$SigmaN
           A <- Sigma[j,j]

           if(lp==0){# no parents
             samplenodescores <- -datatoscore[,j]^2/(2*A) - log(2*pi*A)/2
           } else {
             D <- as.matrix(Sigma[parentnodes,parentnodes])
             choltemp<-chol(D)
             B <- Sigma[j,parentnodes]
             C <- backsolve(choltemp,B,transpose=TRUE)
             E <- backsolve(choltemp,C)

             K <- A - sum(C^2)
             coreMat <- c(1,-E)%*%t(c(1,-E))/K
             xs <- datatoscore[,c(j,parentnodes)]
             samplenodescores <- -rowSums(xs%*%coreMat*xs)/2 - log(2*pi*K)/2
           }
         },
         "bde" = {

           noparams<-2^lp # number of binary states of the parents
           switch(as.character(lp),
                  "0"={# no parents
                    N1<-sum(param$d1[,j],na.rm=TRUE)
                    N0<-sum(param$d0[,j],na.rm=TRUE)
                    NT<-N0+N1
                    theta<-(N1+param$chi/(2*noparams))/(NT+param$chi/noparams) # the probability of each state
                    samplenodescores[which(datatoscore[,j]==1)]<-log(theta) # log scores of 1s
                    samplenodescores[which(datatoscore[,j]==0)]<-log(1-theta) # log scores of 0s
                  },
                  "1"={# one parent

                    summys<-param$data[,parentnodes]
                    summysfull<-datatoscore[,parentnodes]

                    for(i in 1:noparams-1){
                      totest<-which(summys==i)
                      N1<-sum(param$d1[totest,j],na.rm=TRUE)
                      N0<-sum(param$d0[totest,j],na.rm=TRUE)
                      NT<-N0+N1
                      theta<-(N1+param$chi/(2*noparams))/(NT+param$chi/noparams) # the probability of each state
                      toscore<-which(summysfull==i)
                      samplenodescores[toscore[which(datatoscore[toscore,j]==1)]]<-log(theta) # log scores of 1s
                      samplenodescores[toscore[which(datatoscore[toscore,j]==0)]]<-log(1-theta) # log scores of 0s
                    }
                  },
                  { # more parents

                    summys<-colSums(2^(c(0:(lp-1)))*t(param$data[,parentnodes]))
                    tokeep<-which(!is.na(summys+param$d1[,j])) # remove NAs either in the parents or the child
                    if(length(tokeep)<length(summys)){
                      N1s<-collectC(summys[tokeep],param$d1[tokeep,j],noparams)
                      N0s<-collectC(summys[tokeep],param$d0[tokeep,j],noparams)
                    } else {
                      N1s<-collectC(summys,param$d1[,j],noparams)
                      N0s<-collectC(summys,param$d0[,j],noparams)
                    }
                    NTs<-N0s+N1s
                    thetas<-(N1s+param$chi/(2*noparams))/(NTs+param$chi/noparams) # the probability of each state

                    summysfull<-colSums(2^(c(0:(lp-1)))*t(datatoscore[,parentnodes]))
                    ones<-which(datatoscore[,j]==1)
                    samplenodescores[ones]<-log(thetas[summysfull[ones]+1])
                    zeros<-which(datatoscore[,j]==0)
                    samplenodescores[zeros]<-log(1-thetas[summysfull[zeros]+1])
                  })
         },
         "bdecat" = {
           lp<-length(parentnodes) # number of parents
           chi<-param$chi

           Cj <- param$Cvec[j] # number of levels of j

           # Get parameters
           switch(as.character(lp),
                  "0"={# no parents
                    Cp <- 1 # effectively 1 parent level
                    summys <- rep(0, nrow(param$data))
                  },
                  "1"={# one parent
                    Cp <- param$Cvec[parentnodes] # number of parent levels
                    summys <- param$data[,parentnodes]
                  },
                  { # more parents
                    Cp <- prod(param$Cvec[parentnodes])
                    # use mixed radix mapping to unique parent states
                    summys<-colSums(cumprod(c(1,param$Cvec[parentnodes[-lp]]))*t(param$data[,parentnodes]))
                  })

           if(!is.null(param$weightvector)){
             Ns <- BiDAG:::collectCcatwt(summys, param$data[,j], param$weightvector, Cp, Cj)
           } else{
             Ns <- BiDAG:::collectCcat(summys, param$data[,j], Cp, Cj)
           }
           NTs <- rowSums(Ns)

           # Score data
           switch(as.character(lp),
                  "0"={# no parents
                    samplenodescores <- log(Ns[datatoscore[,j]+1] + chi/Cj) - log(NTs + chi)
                  },
                  "1"={# one parent
                    pa_idx <- datatoscore[,parentnodes] + 1
                    j_pa_idx <- Cp * datatoscore[,j] + pa_idx
                    samplenodescores <- log(Ns[j_pa_idx]+chi/(Cp*Cj)) - log(NTs[pa_idx] + chi/Cp)
                  },
                  { # more parents
                    pa_idx <- colSums(cumprod(c(1,param$Cvec[parentnodes[-lp]]))*t(datatoscore[,parentnodes])) + 1
                    j_pa_idx <- Cp * datatoscore[,j] + pa_idx
                    samplenodescores <- log(Ns[j_pa_idx]+chi/(Cp*Cj)) - log(NTs[pa_idx] + chi/Cp)
                  })
         },
         { # pcart case -- approximation
           for (i in c(1:nrow(datatoscore))) {

             data_plus1 <- rbind(param$data, datatoscore[i,])
             score_num <- opt.pcart(data_plus1, parentnodes, j, param$preLevels,
                                    alpha = param$pcart_alpha, kappa = param$pcart_kappa,
                                    response_type = param$response_type)$dataScore
             score_denom <- opt.pcart(param$data, parentnodes, j, param$preLevels,
                                      alpha = param$pcart_alpha, kappa = param$pcart_kappa,
                                      response_type = param$response_type)$dataScore
             samplenodescores[i] <- score_num - score_denom

           }
         })

  return(samplenodescores)
}

# This function scores a data vector against a DAG by marginalising
# over the posterior distribution of the BGe score
# This is equivalent to the difference in scores of the DAG with
# and without the extra data vector

scoreagainstDAGmargBGe <- function(n, scorepar, incidence, datatoscore){

  baselinescore <- DAGscore(scorepar, incidence)

  scorepar2 <- scorepar # store updated scoring components
  scorepar2$N <- scorepar$N + 1 # updated size
  scorepar2$awpN <- scorepar$awpN + 1 # updated parameter

  # update the constant part of the score
  scorepar2$scoreconstvec <- scorepar$scoreconstvec -
    log(pi)/2 + log(scorepar2$am+scorepar2$N-1)/2 - log(scorepar2$am+scorepar2$N)/2 +
    lgamma((1:n-n+scorepar2$awpN)/2) - lgamma((1:n-n-1+scorepar2$awpN)/2)

  # we store part of the score including the T0 matrix and the data covariance matrix
  T0cov <- scorepar$TN - ((scorepar$am*scorepar$N)/(scorepar$am+scorepar$N))* (scorepar$means)%*%t(scorepar$means)

  samplescores <- rep(0, nrow(datatoscore))

  for(ii in 1:nrow(datatoscore)){
    xs <- as.numeric(datatoscore[ii,]) # this has been recentered by subtracting the data means
    scorepar2$means <- scorepar$means + xs/(scorepar$N+1) # update the mean and posterior matrix
    scorepar2$TN <- T0cov + xs%*%t(xs)*scorepar$N/(scorepar$N+1) + (((scorepar2$am)*scorepar2$N)/(scorepar2$am+scorepar2$N)) * (scorepar2$means)%*%t(scorepar2$means)
    samplescores[ii] <- DAGscore(scorepar2, incidence)
  }

  return(samplescores-baselinescore)
}

