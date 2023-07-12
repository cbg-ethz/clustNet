#' @title bestAICsearch
#'
#' @description best AIC search
#'
#' @param binaryMatrix Data to be clustered
#' @param minK Min number of clusters
#' @param maxK Max number of clusters
#' @param chiVec Vector of chi values
#' @param startseed Seed
#' @param nIterations Number of iterations
#' @param AICrange AIC range
#' @param plot_heatmap TRUE if plotting directly
#'
#' @return list of AIC scrores
#' @export
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom grDevices rgb
#'
bestAICsearch <- function(binaryMatrix, minK = 2, maxK = 5, chiVec = c(1e-3,0.5,1,2,3), startseed = 100, nIterations = 50, AICrange = 100, plot_heatmap=TRUE) {
  #
  #     require('ggplot2')
  #     require('reshape2')

  # the following line is to suppress irrelevant notes caused by ggplot
  Var1=Var2=value=aic=0

  ### Check input parameters
  # if(missing(maxK)) stop("Need to input maximum number of clusters k.")
  # if(missing(chiVec)) stop("Need to input a range of chi values.")
  if (!identical(maxK - floor(maxK), 0) || maxK < 2) stop("maxK has to be a positive whole number > 1.")
  if (!identical(minK - floor(minK), 0) || minK < 1) stop("minK has to be a positive whole number > 1.")
  if(any(chiVec<0)) stop("Chi needs to be positive.")

  ### Prepare matrix for plotting
  nrow <- maxK - minK + 1
  ncol <- length(chiVec)
  aics<-matrix(NA,nrow=nrow,ncol=ncol)
  output <- list()

  ### Populate matrix wich AICs
  for(jj in 1:ncol){
    chi<- chiVec[jj]

    print(paste("Currently clustering for chi =", chi))

    if(chi==0){
      chi<-1e-3 # replace 0 by 1e-3 to avoid taking log of 0
    }
    for(kk in minK:maxK){
      bestCluster <- BMMclusterEM(binaryMatrix = binaryMatrix,
                                   chi = chi, k_clust = kk,
                                   startseed = startseed,
                                   nIterations = nIterations)

      if(length(table(bestCluster$newclustermembership))==kk){
        aics[kk - minK + 1,jj]<-bestCluster$testAIC
        output[[kk - minK + 1 + (jj - 1)*(maxK-minK+1)]] <- bestCluster
      }
    }
  }

  plain_aics <- aics

  ks<-rep(0,ncol)
  minaics<-rep(0,ncol)

  minaics <- apply(aics, 2, min, na.rm=TRUE)

  aics<-t(aics)
  aics<-aics - minaics
  topaics<-AICrange
  aics[which(aics>topaics)]<-topaics
  aicsscaled<-t(t(aics)/colSums(aics))
  # heatmap(t(aics),Rowv=NA, Colv=NA,col = rainbow(256))

  divergy<-aics
  # divergy[which(is.na(divergy))]<-maxxy

  rownames(divergy)<-chiVec
  colnames(divergy)<-c(minK:maxK)

  meltdivergy<-melt(divergy)

  if (plot_heatmap==TRUE){
    ggplot(data = meltdivergy, aes(x=Var1, y=Var2, fill=value)) +
       geom_tile()

    middycol<-c(0.8,0.2,0)

    ggheatmap<-ggplot(data = meltdivergy, aes(Var1, Var2, fill = value))+
    # ggheatmap<-ggplot(data = meltdivergy)+
        geom_tile() +
        xlab(expression(chi)) +
        ylab("k") +
        scale_fill_gradient2(high =rgb(0.98,0.98,1), low = "#117777",
                             mid="#88BBBB",space="Lab",na.value="grey75",
                             midpoint=topaics/2,limit = c(0,topaics), name="AIC\nchange\n") +
        scale_y_continuous(breaks=c(minK:maxK)) +
        theme_minimal() +
        theme(axis.title.x = element_text(vjust=-1),axis.title.y = element_text(angle=0,hjust=-0.5,vjust=0.505)) +
        theme(axis.text.x = element_text(angle = 0, vjust = 0.5,size = 20, hjust = 0.6),
              axis.text.y = element_text(angle = 0, vjust = 0.5,size = 20, hjust = 1),
              legend.text=element_text(size=20),
              axis.title=element_text(size=30),
              legend.title=element_text(size=24))+theme(legend.key.size = unit(2,"line")) +
        theme(plot.margin=unit(c(-0.3,-0.3,0.4,0.4),"cm"))

    print(ggheatmap)

    # pdf(paste("heatmapaic.pdf",sep=""), width=7.5, height=6, onefile=F, pointsize=10,  paper="special")

    ggheatmap +
        theme(
            #axis.title.x = element_text("prune and reattach probability"),
            #axis.title.y = element_text("swap two nodes probability"),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.ticks = element_blank())
  }

  # dev.off()
  return(list(output=output, aics=plain_aics, delta_aic=aic))
}


# #' @title bestAICsearch
# #'
# #' @description best AIC search
# #'
# #' @param binaryMatrix Data to be clustered
# #' @param minK Min number of clusters
# #' @param maxK Max number of clusters
# #' @param chiVec Vector of chi values
# #' @param startseed Seed
# #' @param nIterations Number of iterations
# #' @param AICrange AIC range
# #'
# #' @return list of AIC scrores
# #' @export
# #'
# #' @import ggplot2
# #' @importFrom reshape2 melt
# #' @importFrom grDevices rgb
# #'
# bestAICsearch <- function(binaryMatrix, minK = 2, maxK=5, chiVec=c(1e-3,0.5,1,2,3), startseed = 100, nIterations = 50,AICrange=100) {
# #
# #     require('ggplot2')
# #     require('reshape2')
#
#     ### Check input parameters
#     if(missing(maxK)) stop("Need to input maximum number of clusters k.")
#     if(missing(chiVec)) stop("Need to input a range of chi values.")
#     if (!identical(maxK - floor(maxK), 0) || maxK < 2) stop("maxK has to be a positive whole number > 1.")
#     if (!identical(minK - floor(minK), 0) || minK < 1) stop("minK has to be a positive whole number > 1.")
#     if(any(chiVec<0)) stop("Chi needs to be positive.")
#
#     ### Prepare matrix for plotting
#     nrow <- maxK - minK + 1
#     ncol <- length(chiVec)
#     aics<-matrix(NA,nrow=nrow,ncol=ncol)
#     output <- list()
#
#     ### Populate matrix wich AICs
#     for(jj in 1:ncol){
#         chi<- chiVec[jj]
#
#         print(paste("Currently clustering for chi =", chi))
#
#         if(chi==0){
#             chi<-1e-3 # replace 0 by 1e-3 to avoid taking log of 0
#         }
#         for(kk in minK:maxK){
#             bestCluster <- BBMMclusterEM(binaryMatrix = binaryMatrix,
#                                                         chi = chi, k_clust = kk,
#                                                         startseed = startseed,
#                                                         nIterations = nIterations)
#
#             if(length(table(bestCluster$newclustermembership))==kk){
#                 aics[kk - minK + 1,jj]<-bestCluster$testAIC
#                 output[[kk - minK + 1 + (jj - 1)*(maxK-minK+1)]] <- bestCluster
#             }
#         }
#     }
#
#     ks<-rep(0,ncol)
#     minaics<-rep(0,ncol)
#
#     minaics <- apply(aics, 2, min, na.rm=TRUE)
#
#     aics<-t(aics)
#     aics<-aics - minaics
#     topaics<-AICrange
#     aics[which(aics>topaics)]<-topaics
#     aicsscaled<-t(t(aics)/colSums(aics))
#     #heatmap(t(aics),Rowv=NA, Colv=NA,col = rainbow(256))
#
#     divergy<-aics
#     #divergy[which(is.na(divergy))]<-maxxy
#
#     rownames(divergy)<-chiVec
#     colnames(divergy)<-c(minK:maxK)
#
#     meltdivergy<-melt(divergy)
#
#     #ggplot(data = meltdivergy, aes(x=Var1, y=Var2, fill=value)) +
#     #    geom_tile()
#
#     #middycol<-c(0.8,0.2,0)
#
#     # ggheatmap<-ggplot(data = meltdivergy, aes(meltdivergy$Var1, meltdivergy$Var2, fill = meltdivergy$value))+
#     # # ggheatmap<-ggplot(data = meltdivergy)+
#     #     geom_tile() +
#     #     xlab(expression(chi)) +
#     #     ylab("k") +
#     #     scale_fill_gradient2(high =rgb(0.98,0.98,1), low = rgb(0,0.35,0.8),
#     #                          mid=rgb(0.49,0.665,0.9),space="Lab",na.value="grey75",
#     #                          midpoint=topaics/2,limit = c(0,topaics), name="AIC\nchange\n") +
#     #     scale_y_continuous(breaks=c(minK:maxK)) +
#     #     theme_minimal() +
#     #     theme(axis.title.x = element_text(vjust=-1),axis.title.y = element_text(angle=0,hjust=-0.5,vjust=0.505)) +
#     #     theme(axis.text.x = element_text(angle = 0, vjust = 0.5,size = 20, hjust = 0.6),
#     #           axis.text.y = element_text(angle = 0, vjust = 0.5,size = 20, hjust = 1),
#     #           legend.text=element_text(size=20),
#     #           axis.title=element_text(size=30),
#     #           legend.title=element_text(size=24))+theme(legend.key.size = unit(2,"line")) +
#     #     theme(plot.margin=unit(c(-0.3,-0.3,0.4,0.4),"cm"))
#     #
#     # print(ggheatmap)
#     #
#     # #pdf(paste("heatmapaic.pdf",sep=""), width=7.5, height=6, onefile=F, pointsize=10,  paper="special")
#     #
#     # ggheatmap +
#     #     theme(
#     #         #axis.title.x = element_text("prune and reattach probability"),
#     #         #axis.title.y = element_text("swap two nodes probability"),
#     #         panel.grid.major = element_blank(),
#     #         panel.border = element_blank(),
#     #         panel.background = element_blank(),
#     #         axis.ticks = element_blank())
#
#     #dev.off()
#     return(output)
# }

