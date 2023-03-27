graphClust: Network-based clustering with covariate adjustment
-----------

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

`graphClust` is an R package for network-based clustering of binary data using a Bayesian network mixture model and optional covariate adjustment.

<img src="https://github.com/fritzbayer/graphClust/blob/main/vignettes/figures/graphClust3.png" width="1028"/>

Installation
-----------

In order to install the package, it suffices to launch
`R CMD INSTALL path/to/graphClust`
from a terminal, or `make install` from within the package source folder.

Being hosted on GitHub, it is possible to use the `install_github`
tool from an R session:

```{r eval=FALSE}
library("devtools")
install_github("fritzbayer/graphClust", auth_token="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")
```

`graphClust` requires R `>= 3.5`, and depends on 
`BiDAG` (>= 2.0.2), `reshape2`, `pcalg`,
`RBGL`, `clue` and `grDevices`.


Example
-------

```{r eval=FALSE}
library(graphClust)

# Simulate binary data from 3 clusters
k_clust <- 3
ss <- c(400, 500, 600) # samples in each cluster
simulation_data <- sampleData(k_clust = k_clust, n_vars = 20, n_samples = ss)
sampled_data <- simulation_data$sampled_data

# Network-based clustering
cluster_res <- get_clusters(sampled_data, k_clust = k_clust)

# Calculate the AIC 
library(mclust)
adjustedRandIndex(simulation_data$cluster_membership, cluster_res_t$clustermembership)

# Visualize the networks
library(ggplot2)
library(ggraph)
library(igraph)
library(ggpubr)

graphClust::plot_clusters(cluster_res_t)

# Visualize a single network
my_graph <- igraph::graph_from_adjacency_matrix(cluster_res_t$DAGs[[1]], mode="directed")
graphClust::nice_DAG_plot(my_graph)

```
