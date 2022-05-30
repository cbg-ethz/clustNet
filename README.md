netClust: Network-based clustering with causal covariate adjustment
-----------

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

`netClust` is an R package for network-based clustering of binary data using a Bayesian network mixture model and optional covariate adjustment.

<img src="https://github.com/fritzbayer/netClust/blob/main/vignettes/figures/netClust3.png" width="1028"/>

Installation
-----------

In order to install the package, it suffices to launch
`R CMD INSTALL path/to/netClust`
from a terminal, or `make install` from within the package source folder.

Being hosted on GitHub, it is possible to use the `install_github`
tool from an R session:

```{r eval=FALSE}
library("devtools")
install_github("fritzbayer/netClust", auth_token="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")
```

`netClust` requires R `>= 3.5`, and depends on 
`BiDAG` (>= 2.0.2), `reshape2`, `pcalg`,
`RBGL`, `clue` and `grDevices`.


Example
-------

```{r eval=FALSE}
library(netClust)

# Simulate binary data from 3 clusters
ss <- c(400,500,600) # samples in each cluster
myData <- sampleData(kclust = 3, Nvars = 20, sseed = 2, samplesizes = ss)

# Network-based clustering
clusterRes <- netCluster(myData, kclust = 3)
```
