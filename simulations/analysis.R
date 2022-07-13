library(ggpubr)
library(RColorBrewer)
library(netClust)
# library(extrafont)


setwd("/Users/frbayer/Documents/phd_main/packages/netClust/simulations/results")

color_list <- RColorBrewer::brewer.pal(n = 4, name = "RdYlBu")
cols <- color_list[c(1, 3, 4, 3, 4, 3, 4, 3, 4)]

## create legends

# get label
cols <- color_list[c(3, 4, 3, 4, 3, 4, 3, 4)]
res3 <- readRDS("results--k_clust-6--n_vars-20--n_bg-4--n_it-30--different--TRUE.rds")
samples <- res3$correct_samples[,1:2]
test_data <- melt(samples)
names(test_data) <- c("Var1", "Data", "value")

levels(test_data$Data)[1] <- "Variables & Covariates"
levels(test_data$Data)[2] <- "Variables"

my_legend1 <- ggplot2::ggplot(data=test_data, aes(x=reorder(Data, -value), y=value, fill = Data, colour = Data))+
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
  geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
  # geom_point(position = position_jitterdodge(), alpha = 0.3, pch = 21, cex = 0.9) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  labs(y = "Adjusted Rand Index (ARI)", x = "Method") +
  theme_minimal()

legend_1 <- get_legend(my_legend1)

# get label
cols <- color_list[c(1, 3, 4, 3, 4, 3, 4, 3, 4)]
samples <- res3$correct_samples[,1:3]
test_data <- melt(samples)
names(test_data) <- c("Var1", "Data", "value")

levels(test_data$Data)[1] <- "Covariate-Adjustment"
levels(test_data$Data)[2] <- "Variables & Covariates"
levels(test_data$Data)[3] <- "Variables"

my_legend2 <- ggplot2::ggplot(data=test_data, aes(x=reorder(Data, -value), y=value, fill = Data, colour = Data))+
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
  geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
  # geom_point(position = position_jitterdodge(), alpha = 0.3, pch = 21, cex = 0.9) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  labs(y = "Adjusted Rand Index (ARI)", x = "Method") +
  theme_minimal()

legend_2 <- get_legend(my_legend2)



## plot netClust versions for different numbers of covariates
res1 <- readRDS("results--k_clust-6--n_vars-20--n_bg-0--n_it-30--different--TRUE.rds")
nice_plot(res1$correct_samples)
nbg_data <- melt(res1$correct_samples[,1:3])
nbg_data$nbg <- 0
res2 <- readRDS("results--k_clust-6--n_vars-20--n_bg-2--n_it-30--different--TRUE.rds")
nice_plot(res2$correct_samples)
nbg_data_temp <- melt(res2$correct_samples[,1:3])
nbg_data_temp$nbg <- 2
nbg_data <- rbind(nbg_data,nbg_data_temp)
res3 <- readRDS("results--k_clust-6--n_vars-20--n_bg-4--n_it-30--different--TRUE.rds")
nbg_data_temp <- melt(res3$correct_samples[,1:3])
nice_plot1 <- nice_plot(res3$correct_samples)
nbg_data_temp$nbg <- 4
nbg_data <- rbind(nbg_data,nbg_data_temp)
res4 <- readRDS("results--k_clust-6--n_vars-20--n_bg-6--n_it-30--different--TRUE.rds")
nice_plot(res4$correct_samples)
nbg_data_temp <- melt(res4$correct_samples[,1:3])
nbg_data_temp$nbg <- 6
nbg_data <- rbind(nbg_data,nbg_data_temp)
# res5 <- readRDS("results--k_clust-6--n_vars-20--n_bg-8--n_it-30--different--TRUE.rds")
# nice_plot(res5$correct_samples)
# nbg_data_temp <- melt(res5$correct_samples[,1:3])
# nbg_data_temp$nbg <- 8
# nbg_data <- rbind(nbg_data,nbg_data_temp)
# res6 <- readRDS("results--k_clust-6--n_vars-20--n_bg-10--n_it-30--different--TRUE.rds")
# nice_plot(res6$correct_samples)
# nbg_data_temp <- melt(res6$correct_samples[,1:3])
# nbg_data_temp$nbg <- 10
# nbg_data <- rbind(nbg_data,nbg_data_temp)
# res7 <- readRDS("results--k_clust-6--n_vars-20--n_bg-14--n_it-30--different--TRUE.rds")
# nice_plot(res7$correct_samples)
# nbg_data_temp <- melt(res7$correct_samples[,1:3])
# nbg_data_temp$nbg <- 14
# nbg_data <- rbind(nbg_data,nbg_data_temp)



# save
png("~/Desktop/plot_nice.png", width = 14, height = 12, units = 'cm', res = 300)
nice_plot1
#don't forget to embed fonts.
#embed_fonts("plotname.pdf", outfile = "plotname_embed.pdf")
dev.off()


colnames(nbg_data)[2] <- "Method"

plot_nbg <- ggplot2::ggplot(data=nbg_data, aes(x=factor(nbg), y=value, fill = Method, colour = Method))+
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
  # geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
  geom_point(position = position_jitterdodge(), alpha = 0.3, pch = 21, cex = 0.9) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  labs(y = "Adjusted Rand Index (ARI)", x = "Number of Covariates") +
  theme_minimal(); plot_nbg


# save
png("~/Desktop/plot_nbg.png", width = 14, height = 12, units = 'cm', res = 300)
plot_nbg
#don't forget to embed fonts.
#embed_fonts("plotname.pdf", outfile = "plotname_embed.pdf")
dev.off()

# pdf("~/Desktop/plot_nbg.pdf", width = 18, height = 21)
# plot_nbg
# #don't forget to embed fonts.
# #embed_fonts("plotname.pdf", outfile = "plotname_embed.pdf")
# dev.off()
# embed_fonts("~/Desktop/plot_nbg.pdf", outfile="plot_cm_embed.pdf")

# save
png("~/Desktop/plot_graphs.png", width = 18, height = 21, units = 'cm', res = 300)
par(mfrow = c(2, 3))
n_bg <- 4
node_col <- rep("darkgrey", n_bg)
names(node_col) <- as.character(20+1:n_bg)
graph::plot(res3$sampled_results$bayes_nets[[1]]$DAG, nodeAttrs=list(fillcolor=node_col))
graph::plot(res3$sampled_results$bayes_nets[[2]]$DAG, nodeAttrs=list(fillcolor=node_col))
graph::plot(res3$sampled_results$bayes_nets[[3]]$DAG, nodeAttrs=list(fillcolor=node_col))
graph::plot(res3$sampled_results$bayes_nets[[4]]$DAG, nodeAttrs=list(fillcolor=node_col))
graph::plot(res3$sampled_results$bayes_nets[[5]]$DAG, nodeAttrs=list(fillcolor=node_col))
graph::plot(res3$sampled_results$bayes_nets[[6]]$DAG, nodeAttrs=list(fillcolor=node_col))
#don't forget to embed fonts.
#embed_fonts("plotname.pdf", outfile = "plotname_embed.pdf")
dev.off()
par(mfrow = c(1,1))

# save
png("~/Desktop/plot_graphs_blue.png", width = 18, height = 21, units = 'cm', res = 300)
par(mfrow = c(2, 3))
n_bg <- 4
node_col <- rep(cols[3], n_bg)
names(node_col) <- as.character(20+1:n_bg)
graph::plot(res3$sampled_results$bayes_nets[[1]]$DAG, nodeAttrs=list(fillcolor=node_col))
graph::plot(res3$sampled_results$bayes_nets[[2]]$DAG, nodeAttrs=list(fillcolor=node_col))
graph::plot(res3$sampled_results$bayes_nets[[3]]$DAG, nodeAttrs=list(fillcolor=node_col))
graph::plot(res3$sampled_results$bayes_nets[[4]]$DAG, nodeAttrs=list(fillcolor=node_col))
graph::plot(res3$sampled_results$bayes_nets[[5]]$DAG, nodeAttrs=list(fillcolor=node_col))
graph::plot(res3$sampled_results$bayes_nets[[6]]$DAG, nodeAttrs=list(fillcolor=node_col))
#don't forget to embed fonts.
#embed_fonts("plotname.pdf", outfile = "plotname_embed.pdf")
dev.off()
par(mfrow = c(1,1))

## plot netClust versions for different numbers of covariates

res1 <- readRDS("results--k_clust-6--n_vars-20--n_bg-0--n_it-30--different--TRUE.rds")
nice_plot(res1$correct_samples)
nbg_data <- melt(res1$correct_samples[,2:9])
nbg_data$nbg <- 0
res2 <- readRDS("results--k_clust-6--n_vars-20--n_bg-2--n_it-30--different--TRUE.rds")
nice_plot(res2$correct_samples)
nbg_data_temp <- melt(res2$correct_samples[,2:9])
nbg_data_temp$nbg <- 2
nbg_data <- rbind(nbg_data,nbg_data_temp)
res3 <- readRDS("results--k_clust-6--n_vars-20--n_bg-4--n_it-30--different--TRUE.rds")
nbg_data_temp <- melt(res3$correct_samples[,2:9])
nbg_data_temp$nbg <- 4
nbg_data <- rbind(nbg_data,nbg_data_temp)
nice_plot(res3$correct_samples)
res4 <- readRDS("results--k_clust-6--n_vars-20--n_bg-6--n_it-30--different--TRUE.rds")
nice_plot(res4$correct_samples)
nbg_data_temp <- melt(res4$correct_samples[,2:9])
nbg_data_temp$nbg <- 6
nbg_data <- rbind(nbg_data,nbg_data_temp)

colnames(nbg_data)[2] <- "Method"

ggplot2::ggplot(data=nbg_data, aes(x=reorder(nbg, value), y=value, fill = Method, colour = Method))+
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
  # geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
  geom_point(position = position_jitterdodge(), alpha = 0.3, pch = 21, cex = 0.9) +
  scale_colour_manual(values = cols[2:9]) + scale_fill_manual(values = cols[2:9]) +
  labs(y = "Adjusted Rand Index (ARI)", x = "Number of Covariates") +
  theme_minimal()


samples <- res3$correct_samples[,2:9]
test_data <- melt(samples)
test_data$Var3 <- c(rep("netClust", 2*dim(samples)[1]),
                    rep("kmeans", 2*dim(samples)[1]), rep("Mclust", 2*dim(samples)[1]),
                    rep("BMM", 2*dim(samples)[1]))

color_list <- RColorBrewer::brewer.pal(n = 4, name = "RdYlBu")
cols <- color_list[c(3, 4, 3, 4, 3, 4, 3, 4)]

names(test_data) <- c("Var1", "Method", "value", "Var3")

temp_plot <- ggplot2::ggplot(data=test_data, aes(x=reorder(Var3, -value), y=value, fill = Method, colour = Method))+
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
  geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
  # geom_point(position = position_jitterdodge(), alpha = 0.3, pch = 21, cex = 0.9) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  labs(y = "Adjusted Rand Index (ARI)", x = "Method") +
  theme_minimal()+
  theme(legend.position = "none")


p_other <- ggarrange(temp_plot, legend_1, widths=c(0.7,0.3))

# save
png("~/Desktop/plot_other.png", width = 14, height = 12, units = 'cm', res = 300)
p_other
#don't forget to embed fonts.
#embed_fonts("plotname.pdf", outfile = "plotname_embed.pdf")
dev.off()

p1 <- ggplot2::ggplot(data=nbg_data, aes(x=reorder(nbg, value), y=value, fill = Method, colour = Method))+
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
  # geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
  geom_point(position = position_jitterdodge(), alpha = 0.3, pch = 21, cex = 0.9) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  labs(y = "Adjusted Rand Index (ARI)", x = "Number of Covariates") +
  theme_minimal()+
  theme(legend.position = "none")

p_nbg_other <- ggarrange(p1, legend_1, widths=c(0.85,0.15))

# save
png("~/Desktop/plot_nbg_other.png", width = 18, height = 21, units = 'cm', res = 300)
p_nbg_other
#don't forget to embed fonts.
#embed_fonts("plotname.pdf", outfile = "plotname_embed.pdf")
dev.off()





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
    theme_minimal()+
    theme(legend.position = "none")

  # get label
  cols <- color_list[c(1, 3, 4, 3, 4, 3, 4, 3, 4)]
  samples <- samples[,1:3]
  test_data <- melt(samples)
  names(test_data) <- c("Var1", "Data", "value")

  levels(test_data$Data)[1] <- "Covariate-Adjustment"
  levels(test_data$Data)[2] <- "Variables & Covariates"
  levels(test_data$Data)[3] <- "Variables"

  my_legend2 <- ggplot2::ggplot(data=test_data, aes(x=reorder(Data, -value), y=value, fill = Data, colour = Data))+
    geom_boxplot(outlier.shape = NA, alpha = 0.3) +
    # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
    geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
    # geom_point(position = position_jitterdodge(), alpha = 0.3, pch = 21, cex = 0.9) +
    scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
    labs(y = "Adjusted Rand Index (ARI)", x = "Method") +
    theme_minimal()

  legend_2 <- get_legend(my_legend2)

  p_nice <- ggarrange(p1, legend_2, widths=c(0.65,0.35))

  return(p_nice)
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






