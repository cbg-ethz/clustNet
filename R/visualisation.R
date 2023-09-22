
#' @title nice_DAG_plot
#'
#' @description DAG visualization
#'
#' @param my_DAG DAG
#' @param print_direct print DAG if TRUE
#' @param node_size node size vector
#' @param CPDAG if TRUE, then plot CPDAG instead of DAG
#' @param node_colours node colours
#' @param directed TRUE if nodes should be directed
#' @return A plot of the DAG of class c("gg", "ggplot").
#'
#' @export
#'
# #' @importFrom ggraph ggraph geom_edge_arc geom_node_point geom_node_text circle
#' @importFrom igraph V
#' @importFrom methods as
nice_DAG_plot <- function(my_DAG, print_direct=TRUE, node_size=NULL, CPDAG=TRUE, node_colours="#fdae61", directed=TRUE){

  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop(
      "Package \"ggraph\" must be installed to use this function.",
      call. = FALSE
    )
  }

  # set node size if not specified
  if(is.null(node_size)){
    node_size <- 3.5
  }

  # Transform the DAG to a CPDAG for a given variable selection
  if(CPDAG==TRUE){
    temp_var_selection <- 1:dim(my_DAG)[1]
    my_DAG[temp_var_selection,temp_var_selection] <- as(pcalg::dag2cpdag(as(my_DAG[temp_var_selection,temp_var_selection], "graphNEL")), "matrix")
  }

  my_graph <- igraph::graph_from_adjacency_matrix(my_DAG, mode="directed")

  if(directed==TRUE){
    # my_graph <- igraph::graph_from_adjacency_matrix(my_DAG, mode="directed")
    arrowsize <- 2.3
  }else{
    # my_graph <- igraph::graph_from_adjacency_matrix(my_DAG, mode="undirected")
    arrowsize <- 0
  }

  # add labelling
  number_of_bar=length(my_graph)
  id = seq(1, length(my_graph))
  angle= 360 * (id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  hjust <- ifelse(angle > 90 & angle<270, 1, 0)
  angle <- ifelse(angle > 90 & angle<270, angle+180, angle)
  name <- names(igraph::V(my_graph))
  # Type <- as.factor(grp)


  p1 <- ggraph::ggraph(my_graph, layout="circle")+
    ggraph::geom_edge_arc(arrow = ggplot2::arrow(length = ggplot2::unit(arrowsize, 'mm')),
                  start_cap = ggraph::circle(2.3, 'mm'),
                  end_cap = ggraph::circle(2, 'mm'),
                  edge_colour="black", edge_alpha=0.6, edge_width=0.4, ggplot2::aes(circular=TRUE)) +
    ggraph::geom_node_point(size=node_size, color=node_colours, alpha=0.9) +
    ggraph::geom_node_text(ggplot2::aes(label=paste("    ",name,"    "),
                       angle=angle, hjust=hjust), size=2.3, color="black") +
    ggplot2::theme_void() +
    # theme(
    #   # legend.position="none",
    #   plot.margin=unit(c(0,0,0,0), "null"),
    #   panel.spacing=unit(c(0,0,0,0), "null")
    # ) +
    ggplot2::expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2)) +
    ggplot2::coord_fixed()

  if(print_direct){
    print(p1)
  }
  return(p1)

}


# Compute by entropy
entropy <- function(cat.vect){
  # from https://stats.stackexchange.com/questions/221332/variance-of-a-distribution-of-multi-level-categorical-data
  px  = table(cat.vect)/length(cat.vect)
  lpx = log(px, base=2)
  ent = -sum(px*lpx)
  return(ent)
}


#' @title plot_clusters
#'
#' @description Plot clusters
#'
#' @param cluster_results Cluster results
#' @param node_colours node colours
#' @param scale_entropy if true, entropy measure will be used to determine size of the nodes
#' @param directed TRUE if nodes should be directed
#' @return A summary plot of all cluster networks of class c("gg", "ggplot", "ggarrange").
#'
#' @export
#' @examples
#' \donttest{
#' # Simulate data
#' sampled_data <- sampleData(n_vars = 15, n_bg = 0)$sampled_data
#' # learn clusters
#' cluster_results <- get_clusters(sampled_data)
#' # Load additional pacakges to visualize the networks
#' library(ggplot2)
#' library(ggraph)
#' library(igraph)
#' library(ggpubr)
#' # Visualize networks
#' plot_clusters(cluster_results)
#' }
#' @importFrom igraph graph_from_adjacency_matrix
plot_clusters <- function(cluster_results, node_colours="#fdae61", scale_entropy = FALSE, directed=TRUE){

  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop(
      "Package \"ggpubr\" must be installed to use this function.",
      call. = FALSE
    )
  }

  data <- as.data.frame(cluster_results$data)

  node_size <- NULL
  # if entropy should define the size
  if (!is.null(data)){
    k_clust <- length(cluster_results$DAGs)

    # # calculate overall entropy
    # total_entropies <- sapply(data, function(x) entropy(x))
    # # calculate relative entropy and select variables with the lowest entropy
    # all_var_selected <- c()
    # diff_entropies <- matrix(NA, nrow = k_clust, ncol = ncol(data))
    # for (nn in 1:k_clust){
    #   # select data
    #   cluster_data <- data[cluster_results$clustermembership==nn,]
    #   # calculate entropy in clusters
    #   cluster_entropies <- sapply(cluster_data, function(x) entropy(x))
    #   # calculate relative entropy
    #   diff_entropies[nn,] <- cluster_entropies-total_entropies
    #   # select variables with lowest relative entropy
    #   # var_selection <- order(diff_entropies[nn,])[1:n_net]
    #   # all_var_selected <- unique(c(all_var_selected, var_selection))
    #   # top_entropy[[nn]] <- var_selection
    # }
    #
    # # calculate size
    # all_var_selected <- 1:dim(data)[2]
    # node_size_percentage <- 1-(diff_entropies[,all_var_selected]+abs(min(diff_entropies[,all_var_selected])))/(max(diff_entropies[,all_var_selected])+abs(min(diff_entropies[,all_var_selected])))
    # node_size <- 1.5+node_size_percentage*6
    #
    # # calculate overall entropy
    # total_entropies <- sapply(data, function(x) entropy(x))
    # # calculate relative entropy and select variables with the lowest entropy
    # n_net <- 7 # number of selected variables per cluster
    # all_var_selected <- c()
    # diff_entropies <- matrix(NA, nrow = k_clust, ncol = ncol(data))
    # binary_frequency <- matrix(NA, nrow = k_clust, ncol = ncol(data))
    # for (nn in 1:k_clust){
    #   # select data
    #   cluster_data <- data[cluster_results$clustermembership==nn,]
    #   # calculate entropy in clusters
    #   cluster_entropies <- sapply(cluster_data, function(x) entropy(x))
    #   # calculate relative entropy
    #   diff_entropies[nn,] <- cluster_entropies-total_entropies
    #   # calculate frequency
    #   binary_frequency[nn,] <- sapply(cluster_data, function(x) sum(x))
    #   # # # select variables with lowest relative entropy
    #   # # var_selection <- order(diff_entropies[nn,])[1:n_net]
    #   # # all_var_selected <- unique(c(all_var_selected, var_selection))
    #   # # select variables with highest frequency
    #   # var_selection <- order(binary_frequency[nn,], decreasing = TRUE)[1:n_net]
    #   # all_var_selected <- unique(c(all_var_selected, var_selection))
    # }
    #
    # # set node size for entropy
    # node_size_percentage <- 1-(diff_entropies[,]+abs(min(diff_entropies[,])))/(max(diff_entropies[,])+abs(min(diff_entropies[,])))
    # node_size <- 1.5+node_size_percentage*6
    #
    # # if there are binary variables, integrate them
    #   # set cols of genomic information and set frequ to zero elsewhere
    #   binary_frequency[,-binary_cols] <- 0
    #
    #   # set node size for frequency
    #   node_size_percentage_frequ <- (binary_frequency[,]+abs(min(binary_frequency[,])))/(max(binary_frequency[,])+abs(min(binary_frequency[,])))
    #   node_size_frequ <- 1.5+node_size_percentage_frequ*6
    #
    #   # merge
    #   node_size[,binary_cols] <- node_size_frequ[,binary_cols]

    ##
    # top_entropy <- list()
    # calculate overall entropy
    total_entropies <- sapply(data, function(x) entropy(x))
    # calculate relative entropy and select variables with the lowest entropy
    # all_var_selected <- c()
    cluster_dim <- c()
    diff_entropies <- matrix(NA, nrow = k_clust, ncol = ncol(data))
    binary_frequency <- matrix(NA, nrow = k_clust, ncol = ncol(data))
    for (nn in 1:k_clust){
      # select data
      cluster_data <- data[cluster_results$clustermembership==nn,]
      cluster_dim[nn] <- dim(cluster_data)[1]
      # calculate entropy in clusters
      cluster_entropies <- sapply(cluster_data, function(x) entropy(x))
      # calculate relative entropy
      diff_entropies[nn,] <- cluster_entropies-total_entropies
      # calculate frequency
      binary_frequency[nn,] <- sapply(cluster_data, function(x) sum(x))
      # # select variables with lowest relative entropy
      # var_selection <- order(diff_entropies[nn,])[1:n_net]
      # all_var_selected <- unique(c(all_var_selected, var_selection))
      # select variables with highest frequency
      # var_selection <- order(binary_frequency[nn,], decreasing = TRUE)[1:n_net]
      # all_var_selected <- unique(c(all_var_selected, var_selection))
      # top_entropy[[nn]] <- var_selection
    }

    # normalize by cluster size
    binary_frequency <- sweep(binary_frequency, 1, cluster_dim, FUN = "/")

    binary_frequency_temp <- binary_frequency

    # set cols of genomic information and set frequ to zero elsewhere
    # binary_cols <- 1:46
    # binary_frequency[,-binary_cols] <- 0

    # node_size_percentage_frequ <- (binary_frequency[,]+abs(min(binary_frequency[,])))/(max(binary_frequency[,])+abs(min(binary_frequency[,])))
    # node_size_frequ <- 1.5+node_size_percentage_frequ*6
    #
    # node_size[,binary_cols] <- node_size_frequ[,binary_cols]

    # Find the maximum for each column and normalize per row for covariates
    max_col_values <- apply(binary_frequency_temp, 2, max)
    node_size_percentage_frequ_temp <- sweep(binary_frequency_temp, 2, max_col_values, FUN = "/")
    node_size <- 1.5+node_size_percentage_frequ_temp*4

    if(scale_entropy){
      # set node size
      node_size_percentage <- 1-(diff_entropies[,]+abs(min(diff_entropies[,])))/(max(diff_entropies[,])+abs(min(diff_entropies[,])))
      node_size <- 1.5+node_size_percentage*6
    }
  }

  # plot DAGs of each cluster
  p_list <- list()
  k_clust <- length(cluster_results$DAGs)
  for (ii in 1:k_clust){
    # my_graph <- graph_from_adjacency_matrix(cluster_results$DAGs[ii][[1]], mode="directed")
    my_DAG <- cluster_results$DAGs[ii][[1]]
    # my_graph <- igraph::graph_from_adjacency_matrix(cluster_results$DAGs[ii][[1]], mode="directed")
    p_list[[ii]] <- nice_DAG_plot(my_DAG, print_direct=FALSE, node_size=node_size[ii,], node_colours=node_colours, directed=directed)
  }
  ggpubr::ggarrange(plotlist=p_list, labels = paste("Cluster", LETTERS[1:k_clust]))#, ncol = 2, nrow = 2)
}


# #' @title plot_clusters
# #'
# #' @description Plot clusters
# #'
# #' @param cluster_results Cluster results
# #' @param data data
# #' @param node_colours node colours
# #'
# #' @export
# #' @importFrom igraph graph_from_adjacency_matrix
# #' @importFrom ggpubr ggarrange
# plot_clusters <- function(cluster_results, data=NULL, node_colours="#fdae61"){
#
#   node_size <- NULL
#   # if entropy should define the size
#   if (!is.null(data)){
#     k_clust <- length(cluster_results$DAGs)
#
#     # n_net <- 1
#     # top_entropy <- list()
#     # calculate overall entropy
#     total_entropies <- sapply(data, function(x) entropy(x))
#     # calculate relative entropy and select variables with the lowest entropy
#     all_var_selected <- c()
#     diff_entropies <- matrix(NA, nrow = k_clust, ncol = ncol(data))
#     for (nn in 1:k_clust){
#       # select data
#       cluster_data <- data[cluster_results$clustermembership==nn,]
#       # calculate entropy in clusters
#       cluster_entropies <- sapply(cluster_data, function(x) entropy(x))
#       # calculate relative entropy
#       diff_entropies[nn,] <- cluster_entropies-total_entropies
#       # select variables with lowest relative entropy
#       # var_selection <- order(diff_entropies[nn,])[1:n_net]
#       # all_var_selected <- unique(c(all_var_selected, var_selection))
#       # top_entropy[[nn]] <- var_selection
#     }
#
#     # calculate size
#     all_var_selected <- 1:dim(data)[2]
#     node_size_percentage <- 1-(diff_entropies[,all_var_selected]+abs(min(diff_entropies[,all_var_selected])))/(max(diff_entropies[,all_var_selected])+abs(min(diff_entropies[,all_var_selected])))
#     node_size <- 1.5+node_size_percentage*6
#   }
#
#   # plot DAGs of each cluster
#   p_list <- list()
#   k_clust <- length(cluster_results$DAGs)
#   for (ii in 1:k_clust){
#     # my_graph <- graph_from_adjacency_matrix(cluster_results$DAGs[ii][[1]], mode="directed")
#     my_DAG <- cluster_results$DAGs[ii][[1]]
#     # my_graph <- igraph::graph_from_adjacency_matrix(cluster_results$DAGs[ii][[1]], mode="directed")
#     p_list[[ii]] <- nice_DAG_plot(my_DAG, print_direct=FALSE, node_size=node_size[ii,], node_colours=node_colours)
#   }
#   ggarrange(plotlist=p_list, labels = paste("Cluster", LETTERS[1:k_clust]))#, ncol = 2, nrow = 2)
# }



#' @title density_plot
#'
#' @description Create 2d dimensionality reduction of sample fit to Bayesian network clusters
#'
#' @param cluster_results Cluster results from function get_clusters
#' @param var_selection Selected variables to consider, e.g. c(1:5) for first five only
#' @param colourys A vector specifying the colors of each cluster (optional)
#' @return A density plot of class recordedplot.
#'
#' @export
#' @examples
#' \donttest{
#' # Simulate data
#' sampled_data <- sampleData(n_vars = 15, n_samples = c(200,200,200))$sampled_data
#' # Learn clusters
#' cluster_results <- get_clusters(sampled_data)
#' # Load additional pacakges to create a 2d dimensionality reduction
#' library(car)
#' library(ks)
#' library(ggplot2)
#' library(graphics)
#' library(stats)
#' # Plot a 2d dimensionality reduction
#' density_plot(cluster_results)
#' }
density_plot <- function(cluster_results, var_selection = NULL, colourys = NULL){

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




  if (is.null(colourys)){
    colourys<-c("#202020","#774411","#DDAA77","#ed2124","#114477","#CC99BB", "#88CCAA","#117744","#77AADD","#771122","#AA4455","#DD7788","#AA7744",
                         "#777711","#AAAA44","#DDDD77","#44AA77","#117777","#44AAAA","#77CCCC","#4477AA","#771155","#AA4488")
  }

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

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
