netClust: Network-based clustering with optional covariate adjustment
-----------

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

# Simulate binary data of 3 clusters
ss <- c(400,500,600) # samples in each cluster
myData <- sampleData(kclust = 3, Nvars = 20, sseed = 2, samplesizes = ss)

# Network-based clustering
clusterRes <- netCluster(myData, kclust = 3)
```
