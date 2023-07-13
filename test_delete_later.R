function(j,parentnodes,n,param,datatoscore) {
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
                      N1s<-BiDAG:::collectC(summys[tokeep],param$d1[tokeep,j],noparams)
                      N0s<-BiDAG:::collectC(summys[tokeep],param$d0[tokeep,j],noparams)
                    } else {
                      N1s<-BiDAG:::collectC(summys,param$d1[,j],noparams)
                      N0s<-BiDAG:::collectC(summys,param$d0[,j],noparams)
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


