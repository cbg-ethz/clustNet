library(BiDAG)
library(pcalg)
library(RBGL)
library(clue)
library(weights)

library(BBMMclusterEM)
library(graphClust)


ss <- c(400,500,600)
myData <- sampleData(kclust = 3, Nvars = 20, sseed = 2, samplesizes = ss)
# myData <- readRDS("../testData.rds")

parClusterRes <- graphClusterParallel(myData, kclust = 3)

sum(parClusterRes$clustermembership[1:ss[1]]==3)/ss[1]
sum(parClusterRes$clustermembership[(ss[1]+1):(ss[1]+ss[2])]==2)/ss[2]
sum(parClusterRes$clustermembership[(ss[1]+ss[2]+1):(ss[1]+ss[2]+ss[3])]==1)/ss[3]

clusterRes <- graphCluster(myData, kclust = 3)

sum(clusterRes$clustermembership[1:ss[1]]==3)/ss[1]
sum(clusterRes$clustermembership[(ss[1]+1):(ss[1]+ss[2])]==2)/ss[2]
sum(clusterRes$clustermembership[(ss[1]+ss[2]+1):(ss[1]+ss[2]+ss[3])]==1)/ss[3]

binClust <- BBMMclusterEM(binaryMatrix = myData, chi = 1, kclust = 3, startseed = 1, nIterations = 10, verbose=TRUE)

sum(binClust$newclustermembership[1:ss[1]]==3)/ss[1]
sum(binClust$newclustermembership[(ss[1]+1):(ss[1]+ss[2])]==2)/ss[2]
sum(binClust$newclustermembership[(ss[1]+ss[2]+1):(ss[1]+ss[2]+ss[3])]==1)/ss[3]


