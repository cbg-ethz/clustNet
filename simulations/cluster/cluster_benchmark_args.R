#!/usr/bin/env Rscript

library(combinat)
library(netClust)
library(mclust)

# pass parameters

# k_clust n_vars n_bg n_it bgedges equal_cpt_bg
# 4 20 10 20 different TRUE

args = commandArgs(trailingOnly=TRUE)

k_clust <- as.numeric(args[1])
n_vars <- as.numeric(args[2])
n_bg <- as.numeric(args[3])
n_it <- as.numeric(args[4])
# n_samples <- as.numeric(args[5])
bgedges <- as.character(args[5])
equal_cpt_bg <- as.logical(args[6])

# benchmark
results <- benchmark_methods(k_clust = k_clust, n_vars = n_vars, n_bg = n_bg, n_it = n_it, n_samples = NULL,
                             bgedges = bgedges, equal_cpt_bg = equal_cpt_bg)
# save
saveRDS(results, paste0("results/results--k_clust-", k_clust, "--n_vars-", n_vars, "--n_bg-", n_bg, "--n_it-", n_it, ".rds"))

