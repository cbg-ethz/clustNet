
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
#' @param data data
#' @param node_colours node colours
#' @param binary_cols Input can be a vector of the binary variables; if "TRUE", the binary variables are automatically determined; if NULL, then entropy will be used only
#' @param directed TRUE if nodes should be directed
#'
#' @export
#' @importFrom igraph graph_from_adjacency_matrix
# #' @importFrom ggpubr ggarrange
plot_clusters <- function(cluster_results, data=NULL, node_colours="#fdae61", binary_cols = TRUE, directed=TRUE){
  
  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop(
      "Package \"ggpubr\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
  node_size <- NULL
  # if entropy should define the size
  if (!is.null(data)){
    k_clust <- length(cluster_results$DAGs)
    
    # calculate overall entropy
    total_entropies <- sapply(data, function(x) entropy(x))
    # calculate relative entropy and select variables with the lowest entropy
    all_var_selected <- c()
    diff_entropies <- matrix(NA, nrow = k_clust, ncol = ncol(data))
    for (nn in 1:k_clust){
      # select data
      cluster_data <- data[cluster_results$clustermembership==nn,]
      # calculate entropy in clusters
      cluster_entropies <- sapply(cluster_data, function(x) entropy(x))
      # calculate relative entropy
      diff_entropies[nn,] <- cluster_entropies-total_entropies
      # select variables with lowest relative entropy
      # var_selection <- order(diff_entropies[nn,])[1:n_net]
      # all_var_selected <- unique(c(all_var_selected, var_selection))
      # top_entropy[[nn]] <- var_selection
    }
    
    # calculate size
    all_var_selected <- 1:dim(data)[2]
    node_size_percentage <- 1-(diff_entropies[,all_var_selected]+abs(min(diff_entropies[,all_var_selected])))/(max(diff_entropies[,all_var_selected])+abs(min(diff_entropies[,all_var_selected])))
    node_size <- 1.5+node_size_percentage*6
    
    # calculate overall entropy
    total_entropies <- sapply(data, function(x) entropy(x))
    # calculate relative entropy and select variables with the lowest entropy
    n_net <- 7 # number of selected variables per cluster
    all_var_selected <- c()
    diff_entropies <- matrix(NA, nrow = k_clust, ncol = ncol(data))
    binary_frequency <- matrix(NA, nrow = k_clust, ncol = ncol(data))
    for (nn in 1:k_clust){
      # select data
      cluster_data <- data[cluster_results$clustermembership==nn,]
      # calculate entropy in clusters
      cluster_entropies <- sapply(cluster_data, function(x) entropy(x))
      # calculate relative entropy
      diff_entropies[nn,] <- cluster_entropies-total_entropies
      # calculate frequency
      binary_frequency[nn,] <- sapply(cluster_data, function(x) sum(x))
      # # # select variables with lowest relative entropy
      # # var_selection <- order(diff_entropies[nn,])[1:n_net]
      # # all_var_selected <- unique(c(all_var_selected, var_selection))
      # # select variables with highest frequency
      # var_selection <- order(binary_frequency[nn,], decreasing = TRUE)[1:n_net]
      # all_var_selected <- unique(c(all_var_selected, var_selection))
    }
    
    # set node size for entropy
    node_size_percentage <- 1-(diff_entropies[,]+abs(min(diff_entropies[,])))/(max(diff_entropies[,])+abs(min(diff_entropies[,])))
    node_size <- 1.5+node_size_percentage*6
    
    # if there are binary variables, integrate them
    if(!is.null(binary_cols)){
      if(binary_cols==TRUE){
        binary_cols <- as.vector(which(colSums(cluster_data>1)==0))
      }
      # set cols of genomic information and set frequ to zero elsewhere
      binary_frequency[,-binary_cols] <- 0
      
      # set node size for frequency
      node_size_percentage_frequ <- (binary_frequency[,]+abs(min(binary_frequency[,])))/(max(binary_frequency[,])+abs(min(binary_frequency[,])))
      node_size_frequ <- 1.5+node_size_percentage_frequ*6
      
      # merge
      node_size[,binary_cols] <- node_size_frequ[,binary_cols]
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
