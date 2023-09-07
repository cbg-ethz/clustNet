clustNet: Network-based clustering with covariate adjustment
-----------

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

`clustNet` is an R package for network-based clustering of binary data using a Bayesian network mixture model and optional covariate adjustment.

Installation
-----------

The latest stable version of clustNet is available on CRAN and can be installed with

```{r eval=FALSE}
install.packages("clustNet")
```
from within an R session.

Being hosted on GitHub, it is also possible to use the `install_github` tool from an R session to install the latest development version:

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

# Visualize the networks
library(ggplot2)
library(ggraph)
library(igraph)
library(ggpubr)

plot_clusters(cluster_res)

# Visualize a single network
my_graph <- igraph::graph_from_adjacency_matrix(cluster_res$DAGs[[1]], mode="directed")
nice_DAG_plot(my_graph)

```
