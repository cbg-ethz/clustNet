clustNet: Network-based clustering with covariate adjustment
-----------

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

`clustNet` is an R package for network-based clustering of binary data using a Bayesian network mixture model and optional covariate adjustment.

Installation
-----------

In order to install the package, it suffices to launch
`R CMD INSTALL path/to/clustNet`
from a terminal, or `make install` from within the package source folder.

Being hosted on GitHub, it is possible to use the `install_github`
tool from an R session:

```{r eval=FALSE}
library("devtools")
install_github("fritzbayer/clustNet")
```

`clustNet` requires R `>= 3.5`.


Example
-------

```{r eval=FALSE}
library(clustNet)

# Simulate binary data from 3 clusters
k_clust <- 3
ss <- c(400, 500, 600) # samples in each clusterclustNet
simulation_data <- sampleData(k_clust = k_clust, n_vars = 20, n_samples = ss)
sampled_data <- simulation_data$sampled_data

# Network-based clustering
cluster_res <- get_clusters(sampled_data, k_clust = k_clust)

# Calculate the ARI 
library(mclust)
adjustedRandIndex(simulation_data$cluster_membership, cluster_res_t$clustermembership)

# Visualize the networks
library(ggplot2)
library(ggraph)
library(igraph)
library(ggpubr)

clustNet::plot_clusters(cluster_res_t)

# Visualize a single network
my_graph <- igraph::graph_from_adjacency_matrix(cluster_res_t$DAGs[[1]], mode="directed")
clustNet::nice_DAG_plot(my_graph)

```
