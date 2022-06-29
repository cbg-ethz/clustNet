library(combinat)


# simulation of clustering

k_clust <- 5
n_bg <- 5
n_vars <- 10
sseed <- 1

# sample data
sampled_results <- sampleData(k_clust = k_clust, n_vars = n_vars, n_bg = n_bg, sseed = sseed)
sampled_data <- simulation_results$sampled_data

# clustering 



# cluster with covariate-adjusted framework
cluster_results1 <- netClust(sampled_data, kclust = k_clust, nbg = n_bg, EMseeds=1)

correct_samples1 <- max_match(sampled_results$cluster_membership,cluster_results1$clustermembership)

# cluster all variables (variables and covariates)
cluster_results2 <- netClust(sampled_data, kclust = 0, nbg = n_bg, EMseeds=1)

correct_samples2 <- max_match(sampled_results$cluster_membership,cluster_results2$clustermembership)

# cluster only variables without covariates
reduced_data <- sampled_data[,1:n_vars]
cluster_results3 <- netClust(reduced_data, kclust = k_clust, nbg = 0, EMseeds=1)

correct_samples3 <- max_match(sampled_results$cluster_membership,cluster_results3$clustermembership)



correct_fraction <- correct_samples/(dim(sampled_data)[1])



cluster_benchmark <- function(sampled_data, kclust = k_clust, nbg = n_bg, n_vars = n_vars){
  # cluster with covariate-adjusted framework
  cluster_results1 <- netClust(sampled_data, kclust = k_clust, nbg = n_bg, EMseeds=1)
  
  correct_samples1 <- max_match(sampled_results$cluster_membership,cluster_results1$clustermembership)
  
  # cluster all variables (variables and covariates)
  cluster_results2 <- netClust(sampled_data, kclust = 0, nbg = n_bg, EMseeds=1)
  
  correct_samples2 <- max_match(sampled_results$cluster_membership,cluster_results2$clustermembership)
  
  # cluster only variables without covariates
  reduced_data <- sampled_data[,1:n_vars]
  cluster_results3 <- netClust(reduced_data, kclust = k_clust, nbg = 0, EMseeds=1)
  
  correct_samples3 <- max_match(sampled_results$cluster_membership,cluster_results3$clustermembership)
  
  correct_samples <- c(correct_samples1, correct_samples2, correct_samples3)
  correct_fraction <- correct_samples/(dim(sampled_data)[1])
  # correct_fraction <- correct_samples/(dim(sampled_data)[1])
  
}



