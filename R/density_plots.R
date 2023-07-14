
#' @title density_plot
#'
#' @description Create 2d visualisation of sample fit to Bayesian network clusters
#'
#' @param cluster_results Cluster results from function get_clusters
#' @param var_selection Selected variables to consider, e.g. c(1:5) for first five only
#'
#' @export
# #' @importFrom igraph graph_from_adjacency_matrix
density_plot <- function(cluster_results, var_selection = NULL){
  
  if (!requireNamespace("car", quietly = TRUE)) {
    stop(
      "Package \"car\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
  if (!requireNamespace("ks", quietly = TRUE)) {
    stop(
      "Package \"ks\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package \"ggplot2\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
  if (!requireNamespace("graphics", quietly = TRUE)) {
    stop(
      "Package \"graphics\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop(
      "Package \"stats\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
  if (is.null(var_selection)){
    var_selection <- 1:dim(cluster_results$data)[2]
  }
  
  test_data <- cluster_results$data[,var_selection]
  
  # process input data
  clustermembership <- as.factor(cluster_results$clustermembership)
  levels(clustermembership) <- LETTERS[1:length(levels(clustermembership))]
  levelly <- levels(clustermembership)
  
  # cancer_type <- as.factor(cancer_type)
  # levels(cancer_type) <- c("AML", "CMML", "MDS", "MPN")
  # levelly2 <- levels(cancer_type)
  # colourys2 <- c("#DD7788", "#771122", "#117777", "#DDDD77")
  
  k_clust <- length(cluster_results$DAGs)
  
  # Calculate scores against clusters
  scoresagainstclusters <- matrix(NA,dim(test_data)[1],k_clust)
  for (k in 1:k_clust){
    allrelativeprobabs <- rep(0, dim(test_data)[1])
    allrelativeprobabs[cluster_results$clustermembership==k] <- 1
  
    scorepar <- BiDAG::scoreparameters("bdecat", as.data.frame(test_data), weightvector = allrelativeprobabs) #, bdepar = bdepar)
    scoresagainstclusters[,k] <- BiDAG::scoreagainstDAG(scorepar,cluster_results$DAGs[[k]][var_selection,var_selection],test_data, bdecatCvec = apply(test_data, 2, function(x) length(unique(x))))
  }
  
  # Calculate divergence matrix
  matsize <- nrow(test_data)
  divergy <- matrix(0, matsize, matsize)
  for(k1 in 1:(matsize-1)){
    P <- exp(scoresagainstclusters[k1,] - max(scoresagainstclusters[k1,]))
    P <- P/sum(P)
    for(k2 in (k1+1):matsize){
      Q <- exp(scoresagainstclusters[k2,] - max(scoresagainstclusters[k2,]))
      Q <- Q/sum(Q)
      M <- (P+Q)/2
      JS <- (P%*%log(P/M) + Q%*%log(Q/M))/2
      divergy[k1,k2] <- JS
      divergy[k2,k1] <- JS
    }
  }
  
  # Perform multidimensional scaling
  reducey <- stats::cmdscale(divergy, k=2)
  
  # plot(reducey, col=cancer_type+1)
  
  # Set colors and plotting parameters
  # colourys<-c("#202020","#774411","#DDAA77","#ed2124","#114477","#CC99BB",
  #                      "#88CCAA","#117744","#77AADD")
  
  colourys<-c("#202020","#771122","#AA4455","#DD7788","#774411","#AA7744",
                       "#DDAA77","#777711","#AAAA44","#DDDD77","#117744","#44AA77",
                       "#88CCAA","#117777","#44AAAA","#77CCCC","#114477","#4477AA",
                       "#77AADD","#771155","#AA4488","#CC99BB")
  
  
  colourys3<-ggplot2::alpha(colourys,0.3)
  graphics::par(mar = c(0,0,0,0))
  
  xlims<-c(1.08*min(reducey[,1]), max(reducey[,1]))
  ylims<-c(min(reducey[,2]), max(reducey[,2]))
  
  graphics::par(bty = 'n')
  plot(1, type="n", axes=F, xlab="", ylab="",xlim=xlims,ylim=ylims, bty="n")
  
  textheights<-(ylims[2]-ylims[1])*rev(1*((c(1:16)-1.5)/(16-1))) +ylims[1]
  textxs<-xlims[1]
  
  k_clust<-length(levelly)
  
  for(k in 1:k_clust){
    selecteddots<-which(clustermembership==levelly[k])
    combdata<-reducey[selecteddots,]
  
    z <- ks::kde(combdata)
  
    plot(z,display="filled.contour2",add=TRUE,cont=c(50),col=c("transparent",paste0(colourys[k],"66")), alpha=0.3, drawpoints=FALSE,drawlabels=FALSE,lwd=1.5, bty="n")
  
    graphics::points(textxs,textheights[k],col=colourys[k],pch=19,cex=1.5)
    graphics::text(textxs,textheights[k],levelly[k],pos=4)
  
  }
  
  # # select mutations
  # mutation_cols <- 1:46
  # mut_data <- test_data[,mutation_cols]
  #
  # #TP53 check
  # tp53_index <- which(colnames(mut_data)=="TP53")
  # TP53only <- reducey[which((mut_data[,tp53_index]==1)&as.vector((rowSums(mut_data)==1)))[1],]
  # # text(TP53only[1], TP53only[2], "TP53 only")
  #
  # #no muts check
  # nomuts <- reducey[which(as.vector(rowSums(mut_data)==0))[1],]
  #
  # # plot edges for TP53 and no mutations
  # arrows(TP53only[1] + 0.05, TP53only[2] - 0.05, TP53only[1] + 0.005, TP53only[2] - 0.005, length=0.1)
  # text(TP53only[1] + 0.06, TP53only[2] - 0.085, paste0("TP53 only\n(",sum((mut_data[,1]==1)&as.vector((rowSums(mut_data)==1)))," samples)"))
  # arrows(nomuts[1] + 0.05, nomuts[2] + 0.05, nomuts[1] + 0.005, nomuts[2] + 0.005, length=0.1)
  # # text(nomuts[1] - 0.06, nomuts[2] - 0.085, paste0("No mutations\n(",sum(as.vector(rowSums(mut_data)==0))," samples)"))
  # text(nomuts[1] + 0.06, nomuts[2] + 0.085, paste0("No mutations and\n no cytogenetic abnormalities \n(",sum(as.vector(rowSums(mut_data)==0))," samples)"))
  #
  # transparent_plot_plain <- recordPlot()
  
  textxs<-xlims[2]*0.88
  
  k_clust<-length(levelly)
  
  for(k in 1:k_clust){
    selecteddots<-which(clustermembership==levelly[k])
  
    graphics::points(reducey[selecteddots,1],reducey[selecteddots,2],col=colourys[k],pch=19,cex=0.5)
  
    # points(textxs,textheights[k],col=colourys[k],pch=19,cex=1.5)
    # graphics::text(textxs,textheights[k],levelly[k],pos=4)
  
  }
  
  transparent_plot <- grDevices::recordPlot()
  
  return(transparent_plot)
}
