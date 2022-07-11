

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





