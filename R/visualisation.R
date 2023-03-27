
#' @title nice_DAG_plot
#'
#' @description DAG visualization
#'
#' @param my_graph DAG
#' @param print_direct print DAG if TRUE
#'
#' @export
nice_DAG_plot <- function(my_graph, print_direct=TRUE){

  # add labelling
  number_of_bar=length(my_graph)
  id = seq(1, length(my_graph))
  angle= 360 * (id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  hjust <- ifelse(angle > 90 & angle<270, 1, 0)
  angle <- ifelse(angle > 90 & angle<270, angle+180, angle)
  name <- names(V(my_graph))
  # Type <- as.factor(grp)

  p1 <- ggraph(my_graph, layout="circle")+
    geom_edge_arc(arrow = arrow(length = unit(2.3, 'mm')),
                  start_cap = circle(2.3, 'mm'),
                  end_cap = circle(2, 'mm'),
                  edge_colour="black", edge_alpha=0.6, edge_width=0.4, aes(circular=TRUE)) +
    geom_node_point(size=3.5, color="#fdae61", alpha=0.9) +
    geom_node_text(aes(label=paste("    ",name,"    "),
                       angle=angle, hjust=hjust), size=2.3, color="black") +
    theme_void() +
    # theme(
    #   # legend.position="none",
    #   plot.margin=unit(c(0,0,0,0), "null"),
    #   panel.spacing=unit(c(0,0,0,0), "null")
    # ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2)) +
    coord_fixed()

  if(print_direct){
    print(p1)
  }
  return(p1)
}


#' @title plot_clusters
#'
#' @description Plot clusters
#'
#' @param cluster_results Cluster results
#'
#' @export
plot_clusters <- function(cluster_results){
  # plot DAGs of each cluster
  p_list <- list()
  k_clust <- length(cluster_results$DAGs)
  for (ii in 1:k_clust){
    my_graph <- graph_from_adjacency_matrix(cluster_results$DAGs[ii][[1]], mode="directed")
    p_list[[ii]] <- nice_DAG_plot(my_graph, print_direct=FALSE)
  }
  ggarrange(plotlist=p_list, labels = paste("Cluster", LETTERS[1:k_clust]))#, ncol = 2, nrow = 2)
}
