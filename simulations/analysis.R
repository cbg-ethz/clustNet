

setwd("/Users/frbayer/Documents/phd_main/packages/netClust/simulations/results")


res1 <- readRDS("results--k_clust-2--n_vars-20--n_bg-10--n_it-20--different--TRUE.rds")
res4 <- readRDS("results--k_clust-4--n_vars-20--n_bg-10--n_it-20--different--TRUE.rds")
res7 <- readRDS("results--k_clust-6--n_vars-20--n_bg-10--n_it-20--different--TRUE.rds")
res8 <- readRDS("results--k_clust-8--n_vars-20--n_bg-10--n_it-20--different--TRUE.rds")
boxplot(res1$correct_samples);
boxplot(res4$correct_samples);
boxplot(res7$correct_samples);
boxplot(res8$correct_samples)

res2 <- readRDS("results--k_clust-4--n_vars-20--n_bg-5--n_it-30--different--TRUE.rds")
boxplot(res2$correct_samples)

res3 <- readRDS("results--k_clust-4--n_vars-20--n_bg-10--n_it-20--different--FALSE.rds")
res5 <- readRDS("results--k_clust-4--n_vars-20--n_bg-10--n_it-20--same--FALSE.rds")
res6 <- readRDS("results--k_clust-4--n_vars-20--n_bg-10--n_it-20--same--TRUE.rds")
boxplot(res3$correct_samples)
boxplot(res5$correct_samples)
boxplot(res6$correct_samples)

res9 <- readRDS("results--k_clust-4--n_vars-20--n_bg-5--n_it-20--different--TRUE.rds")
res10 <- readRDS("results--k_clust-6--n_vars-20--n_bg-5--n_it-20--different--TRUE.rds")
# res11 <- readRDS("results--k_clust-8--n_vars-20--n_bg-5--n_it-20--different--FALSE.rds")
boxplot(res9$correct_samples)
boxplot(res10$correct_samples)




res2 <- readRDS("results--k_clust-6--n_vars-20--n_bg-2--n_it-30--different--TRUE.rds")
nice_plot(res2$correct_samples)
nbg_data <- melt(res2$correct_samples[,1:3])
nbg_data$nbg <- 2
res3 <- readRDS("results--k_clust-6--n_vars-20--n_bg-4--n_it-30--different--TRUE.rds")
nbg_data_temp <- melt(res3$correct_samples[,1:3])
nbg_data_temp$nbg <- 4
nbg_data <- rbind(nbg_data,nbg_data_temp)
nice_plot(res3$correct_samples)
res4 <- readRDS("results--k_clust-6--n_vars-20--n_bg-6--n_it-30--different--TRUE.rds")
nice_plot(res4$correct_samples)
nbg_data_temp <- melt(res4$correct_samples[,1:3])
nbg_data_temp$nbg <- 6
nbg_data <- rbind(nbg_data,nbg_data_temp)

ggplot2::ggplot(data=nbg_data, aes(x=reorder(nbg, value), y=value, fill = Var2, colour = Var2))+
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
  # geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
  geom_point(position = position_jitterdodge(), alpha = 0.3, pch = 21, cex = 0.9) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  labs(y = "Adjusted Rand Index (ARI)", x = "Number of Covariates") +
  theme_minimal()

nice_plot <- function(samples){
  test_data <- melt(samples)
  test_data$Var3 <- c(rep("Cov-Adjust", dim(samples)[1]), rep("netClust", 2*dim(samples)[1]),
                      rep("kmeans", 2*dim(samples)[1]), rep("Mclust", 2*dim(samples)[1]),
                      rep("BMM", 2*dim(samples)[1]))

  color_list <- RColorBrewer::brewer.pal(n = 4, name = "RdYlBu")
  cols <- color_list[c(1, 3, 4, 3, 4, 3, 4, 3, 4)]

  names(test_data) <- c("Var1", "Method", "value", "Var3")

  p1 <- ggplot2::ggplot(data=test_data, aes(x=reorder(Var3, -value), y=value, fill = Method, colour = Method))+
    geom_boxplot(outlier.shape = NA, alpha = 0.3) +
    # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
    geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
    # geom_point(position = position_jitterdodge(), alpha = 0.3, pch = 21, cex = 0.9) +
    scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
    labs(y = "Adjusted Rand Index (ARI)", x = "Method") +
    theme_minimal()

  return(p1)
}

ggplot2::ggplot(data=test_data, aes(x=reorder(Var3, -value), y=value, fill = Method, colour = Method))+
  geom_boxplot(aes(x=reorder(Var3, -value)), outlier.shape = NA, alpha = 0.3) +
  # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
  geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  labs(y = "NRMSE", x = "Dimensions") +
  theme_minimal()



ggplot2::ggplot(data=test_data, aes(x=reorder(Var3, -value), y=value, fill = Var2, colour = Var2))+
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
  geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  labs(y = "Adjusted Rand Index (ARI)", x = "Method") +
  theme_minimal()

scale_fill_discrete(name="Title", labels=c("1","2","3"))








test_data <- melt(res9$correct_samples)

color_list <- brewer.pal(n = 4, name = "RdYlBu")
cols <- color_list[c(1, 3, 4, 3, 4, 3, 4, 3, 4)]

ggplot(data=test_data, aes(x=Var2, y=value, fill = Var2, colour = Var2))+
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
  geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  labs(y = "NRMSE", x = "Dimensions") +
  theme_minimal()

p1 <- ggplot(data = resAllDims2, aes(x = factor(dimBN,
                                                level = c(labelBP[1], labelBP[2], labelBP[3], labelBP[4])),
                                     y = BoxPlotResNRMS, fill = Method, colour = Method),
             alpha = 0.3) + geom_boxplot(outlier.shape = NA,
                                         alpha = 0.3) + geom_point(pch = 21, data = resAllDims2,
                                                                   position = position_jitterdodge(0.15), cex = 0.4,
                                                                   alpha = 0.3) + scale_colour_manual(values = color_list[c(1,
                                                                                                                            3, 4)]) + scale_fill_manual(values = color_list[c(1,
                                                                                                                                                                              3, 4)]) + labs(y = "NRMSE", x = "Dimensions") +
  ylim(0, 1.5) + theme_minimal() + theme(legend.position = "bottom")
}





