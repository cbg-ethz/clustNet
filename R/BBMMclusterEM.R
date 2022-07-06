BBMMclusterEM <- function(binaryMatrix, chi, kclust, startseed=100, nIterations=50, verbose=FALSE) {

    # set.seed(startseed)
    # Check for input arguments
    if (missing(binaryMatrix) || !all(binaryMatrix < 2)) stop("Need a binary matrix as input to cluster.")
    if(missing(chi)) stop('Need to provide a value for chi.')
    if(chi==0) {
        #print('Zero chi, setting chi to 1e-3')
        chi<-1e-3
    }
    if (missing(kclust)) stop('Need to provide a value for kclust.')
    if (as.integer(nIterations) < 1) {
        stop('Need to specify a positive integer.')
    } else {
        nIterations <- as.integer(nIterations)
    }
    if (startseed < 0) stop("Need to specify a positive integer as startseed.")

    tmp <- lapply(seq(nIterations), doIterate, startseed, chi, kclust, binaryMatrix,verbose)

    idx <- which.min(sapply(tmp, '[', 'testAIC'))

    output <- tmp[[idx]]

    return(output)
}

doIterate <- function(idx, startseed, chi, kclust, datatocluster,verbose) {

    seednumber <- startseed + idx
    # set.seed(seednumber)

    if(verbose==TRUE){
      print(paste("Seed", seednumber, "with", kclust, "clusters and", chi, "pseudocounts"))
    }

    clusterresults <- BBMMclusterEMcore(kclust,chi,datatocluster)

    # finally assign to the maximal cluster
    newclustermembership <- reassignsamples(clusterresults$relativeweights)
    # could also output the relative probability of being in that cluster?

    clustersizes <- table(newclustermembership)
    kfound <- length(which(clustersizes>0))

   # these values include the pseudocounts

    #totalloglike <- calcloglike(clusterresults$scoresagainstclusters,clusterresults$tauvec)

    #totalAIC <- 2*kfound*ncol(datatocluster)+2*(kfound-1)-2*totalloglike
    #totalBIC <- log(nrow(datatocluster))*(kfound*ncol(datatocluster)+kfound-1)-2*totalloglike

    #print(totalloglike)
    #print(totalAIC)
    #print(totalBIC)

    # now we recompute with (almost) no pseudocounts to get a ML limit

    againstclusterresults <- scoreagainstemptyEMcluster(kclust,1e-3,clusterresults$relativeweights,clusterresults$tauvec,datatocluster,datatocluster)

    testloglike <- calcloglike(againstclusterresults$scoresagainstclusters,clusterresults$tauvec)
    #print(testloglike)

    testAIC <- 2*kfound*ncol(datatocluster)+2*(kfound-1)-2*testloglike
    testBIC <- log(nrow(datatocluster))*(kfound*ncol(datatocluster)+kfound-1)-2*testloglike

    #print(testAIC)
    #print(testBIC)

    if(verbose==TRUE){
      print(paste("Log likelihood is",testloglike))
      print(paste("AIC is",testAIC))
      print(paste("BIC is",testBIC))
    }

    return(list(seed = seednumber, testAIC = testAIC, newclustermembership = newclustermembership, relativeweights=clusterresults$relativeweights))
}

# This function samples an element from a vector properly

propersample <- function(x){if(length(x)==1) x else sample(x,1)}

# this function assigns sample to the cluster with the highest weight

reassignsamples <- function(sampleprobs){
    newclustermembership <-apply(sampleprobs,1,propersample(which.max)) # find the max per row
    return(newclustermembership)
}

# this function takes in log scores and returns normalised probabilities

allrelativeprobs <- function(samplescores){
    maxscorey<-apply(samplescores,1,max) # find the max of each row
    relativeprobs<-exp(samplescores-maxscorey) # remove max for numerical stability and exponentiate
    relativeprobs<-relativeprobs/rowSums(relativeprobs) # normalise
    return(relativeprobs)
}

# this function takes in probabilities, weights them by the vector tau
# and returns normalised probabilities

allrelativeweights <- function(sampleprobs,tau){
    relativeprobs<-t(t(sampleprobs)*tau)
    relativeprobs<-relativeprobs/rowSums(relativeprobs) # normalise
    return(relativeprobs)
}

# this function takes in log scores with the weight vector and returns the log likelihood

calcloglike <- function(samplescores,tau){
    maxscorey<-apply(samplescores,1,max) # find the max of each row
    loglike<-sum(log(colSums(t(exp(samplescores-maxscorey))*tau))+maxscorey) # remove max for numerical stability and exponentiate
    return(loglike)
}

BBMMclusterEMcore <- function(kclust, chi, datatocluster){

    nbig<-ncol(datatocluster)
    mbig<-nrow(datatocluster)
    scoresagainstclusters<-matrix(0,mbig,kclust)

    diffystart<-10
    diffy<-diffystart # to check for convergence
    county<-0 # to check how many loops we run
    countlimit<-1e4 # hard limit on the number of loops
    errortol<-1e-10# when to stop the assignment

    clustermembership<-sample.int(kclust,mbig,replace=TRUE) # start with random groupings

    for(k in 1:kclust){
        clustersamps<-which(clustermembership==k) # find the members of the cluster
        scoresagainstclusters[clustersamps,k]<-1e-3 # increase the probabilty of clustermembership
        # to get non-uniform starting point
    }

    # this is the weights of each sample for each cluster
    relativeprobabs<-allrelativeprobs(scoresagainstclusters)
    relativeweights<-relativeprobabs # for the starting value

    # main loop
    while((diffy>0) && (county<countlimit)){

        county<-county+1 # update counter

        # first given the current weights we can update tau

        rowtots<-colSums(relativeweights) + chi # add prior to clustersizes
        tauvec<-rowtots/sum(rowtots)

        # and the posterior means

        for(kk in 1:kclust){
            weightvec<-relativeweights[,kk]
            Datawone<- datatocluster*weightvec # the weighted 1s
            #Datawzero<- (t(1-Datafull)*weightvec) # the weighted 0s
            thetas<-(colSums(Datawone)+0.5*chi)/(sum(weightvec)+chi) # we add the prior
            # the posterior means give the (log) probability of the observations
            scoresagainstclusters[,kk]<-colSums(t(datatocluster)*log(thetas)+t(1-datatocluster)*log(1-thetas))
        }

        # now we can update the weights

        relativeprobabs<-allrelativeprobs(scoresagainstclusters)
        newrelativeweights<-allrelativeweights(relativeprobabs,tauvec)

        # calculate the difference
        diffy3<-sum((newrelativeweights-relativeweights)^2)

        if(diffy3<errortol){
            diffy<-diffy-1
        } else {
            diffy<-diffystart # otherwise reset the counter
        }

        relativeweights<-newrelativeweights # for the next loop
    }

    output<-vector("list",0)
    output$relativeweights<-relativeweights
    output$relativeprobs<-relativeprobabs
    output$scoresagainstclusters<-scoresagainstclusters
    output$tauvec<-tauvec

    return(output)

}

scoreagainstemptyEMcluster <- function(kclust, chi, relativeweights, tauvec, datatocluster, datatoscore){

    nbig<-ncol(datatoscore)
    mbig<-nrow(datatoscore)

    scoresagainstclusters<-matrix(0,mbig,kclust)

    for(kk in 1:kclust){
        weightvec<-relativeweights[,kk]
        Datawone<- datatocluster*weightvec # the weighted 1s
        #Datawzero<- (t(1-Datafull)*weightvec) # the weighted 0s
        thetas<-(colSums(Datawone)+0.5*chi)/(sum(weightvec)+chi) # we add the prior
        # the posterior means give the (log) probability of the observations
        scoresagainstclusters[,kk]<-colSums(t(datatoscore)*log(thetas)+t(1-datatoscore)*log(1-thetas))
    }

    # calculate the weights

    relativeprobabs<-allrelativeprobs(scoresagainstclusters)
    relativeweights<-allrelativeweights(relativeprobabs,tauvec)

    output<-vector("list",0)
    output$relativeweights<-relativeweights
    output$relativeprobs<-relativeprobabs
    output$scoresagainstclusters<-scoresagainstclusters

    return(output)

}
