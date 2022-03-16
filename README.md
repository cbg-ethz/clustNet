graphClust: Graph-based clustering of binary data with optional covariate adjustment
-----------

`graphClust` is an R package for graph-based clustering of binary data using a Bayesian network mixture model and optional covariate adjustment.

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
`BiDAG` (>= 2.0.2), `ggplot2`, `reshape2`, `pcalg`,
`RBGL`, `clue` and `grDevices`.
