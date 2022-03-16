netClust: Network-based clustering with optional covariate adjustment
-----------

`netClust` is an R package for network-based clustering of binary data using a Bayesian network mixture model and optional covariate adjustment.

<img src="https://github.com/fritzbayer/netClust/blob/main/blob/netClust3.png" width="528"/>

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
`BiDAG` (>= 2.0.2), `ggplot2`, `reshape2`, `pcalg`,
`RBGL`, `clue` and `grDevices`.
