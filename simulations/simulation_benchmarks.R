library(combinat)
library(netClust)
library(mclust)

# simulation of clustering

k_clust <- 5
n_bg <- 10
n_vars <- 20
sseed <- 1

# sample data
sampled_results <- netClust:::sampleData(k_clust = k_clust, n_vars = n_vars, n_bg = n_bg, sseed = sseed)
sampled_data <- sampled_results$sampled_data
sampled_membership <- sampled_results$cluster_membership

# clustering

correct_samples <- netClust:::cluster_benchmark(sampled_data, sampled_membership, kclust = k_clust, nbg = n_bg, n_vars = n_vars, n_rep = 2)

# correct_fraction <- correct_samples/(dim(sampled_data)[1])

saveRDS(correct_samples, "benchmark_res.rds")



# k-means
res_kmeans1 <- kmeans(sampled_data, 5)
max_match(sampled_membership, res_kmeans1$cluster)

res_kmeans2 <- kmeans(reduced_data, 5)
max_match(sampled_membership, res_kmeans2$cluster)

# Mclust
res_mclust1 <- Mclust(sampled_data, 5)
max_match(sampled_membership, res_mclust1$classification)

res_mclust2 <- Mclust(reduced_data, 5)
max_match(sampled_membership, res_mclust2$classification)

# Bernoulli Mixture Model (BBMMclusterEM)
res_BBMM1 <- BBMMclusterEM(sampled_data, chi = 0.5, kclust = 5, startseed = 1, nIterations = 10, verbose=TRUE)
max_match(sampled_membership, res_BBMM1$newclustermembership)

res_BBMM2 <- BBMMclusterEM(reduced_data, chi = 0.5, kclust = 5, startseed = 1, nIterations = 10, verbose=TRUE)
max_match(sampled_membership, res_BBMM2$newclustermembership)






# cluster with covariate-adjusted framework
cluster_results1 <- netClust(sampled_data, kclust = k_clust, nbg = n_bg, EMseeds=1)

correct_samples1 <- max_match(sampled_results$cluster_membership,cluster_results1$clustermembership)

# cluster all variables (variables and covariates)
cluster_results2 <- netClust(sampled_data, kclust = k_clust, nbg = 0, EMseeds=1)

correct_samples2 <- max_match(sampled_results$cluster_membership,cluster_results2$clustermembership)

# cluster only variables without covariates
reduced_data <- sampled_data[,1:n_vars]
cluster_results3 <- netClust(reduced_data, kclust = k_clust, nbg = 0, EMseeds=1)

correct_samples3 <- max_match(sampled_results$cluster_membership,cluster_results3$clustermembership)





# correct_fraction <- correct_samples/(dim(sampled_data)[1])




