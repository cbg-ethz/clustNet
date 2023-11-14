clustNet: Network-based clustering with covariate adjustment
-----------

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

`clustNet` is an R package for network-based clustering of categorical data using a Bayesian network mixture model and optional covariate adjustment.

Installation
-----------

The latest stable version of clustNet is available on CRAN and can be installed with

```{r eval=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("Rgraphviz", "RBGL"))

install.packages("clustNet")
```
from within an R session. On a normal computer, this should take around 5-60 seconds, depending on how many of the required packages are already installed.

BiocManager::install("remotes")

Being hosted on GitHub, it is also possible to use the `install_github` tool from an R session to install the latest development version:

```{r eval=FALSE}
library("devtools")
install_github("cbg-ethz/clustNet")
```

`clustNet` requires R `>= 3.5`.


Example
-------

```{r eval=FALSE}
library(clustNet)

# Simulate data
k_clust <- 3 # numer of clusters
ss <- c(400, 500, 600) # samples in each cluster
simulation_data <- sampleData(k_clust = k_clust, n_vars = 20, n_samples = ss)
sampled_data <- simulation_data$sampled_data

# Network-based clustering
cluster_results <- get_clusters(sampled_data, k_clust = k_clust)

# Load additional pacakges to visualize the networks
library(ggplot2)
library(ggraph)
library(igraph)
library(ggpubr)

# Visualize networks
plot_clusters(cluster_results)

# Load additional pacakges to create a 2d dimensionality reduction
library(car)
library(ks)
library(graphics)
library(stats)

# Plot a 2d dimensionality reduction
density_plot(cluster_results)

```

On a normal computer, the clustering should take around 2-4 minutes.
