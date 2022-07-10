

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

