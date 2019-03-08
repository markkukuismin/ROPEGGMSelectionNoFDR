# ROPEGGMSelectionNoFDR

# You need following R packages to run the script:

install.packages("ggplot2")<br/>
install.packages("reshape")<br/>
install.packages("camel")<br/>
install.packages("Matrix")<br/>
install.packages("spatstat")<br/>
install.packages("huge")<br/>
install.packages("igraph")<br/>
install.packages("BDgraph")<br/>
install.packages("GGMridge")<br/>
install.packages("rags2ridges")<br/>
install.packages("FastGGM")<br/>
install.packages("fdrtool")

# Example

Explore network clusters from penalized precision matrix estimates:

```r
library(huge)
library(igraph)
library(rags2ridges)
library(fdrtool)

source("ROPE.txt")
source("crossvalidationROPE.txt")
source("ROPEEEDKernelMethod.r")
source("ROPEEEDEmpiricPoisson.r")
source("Scale.r")

lambda = seq(0.01,10,length.out=40) # The largest value should be even larger when p increases

alpha = 0.05 # 100*(1-alpha)% confidence level

B = NULL # Number of breaks in the histogram. Default "Sturges"

p = 200
n = 100

nmb = p*(p-1)/2

I = diag(1,p)

#Simulate data:
set.seed(242377)
L = huge.generator(n=n,d=p,g=4,graph="cluster")  # Generate data with hub structures

Sigma = L$sigma
Theta = L$theta
Theta = as.matrix(Theta)

TrueA = L$theta
TrueA = as.matrix(TrueA)
TrueG = graph.adjacency(TrueA,mode="undirected",diag=F)

ClustersTrue = walktrap.community(TrueG)
TrueNmbOfClusters = length(table(membership(ClustersTrue)))

X = mvrnorm(2*n,rep(0,p),Sigma)
X = huge.npn(X)

Yvalid = X[1:n,]
Y = X[(n+1):(2*n),]

# Choose optimal value for the tuning parameter from the validation set

lambdaROPE = ROPEPoissonEED(Yvalid,lambda=lambda,alpha=alpha,target=I)$lambda

ROPEPoisson = ROPEPoissonEED(Y,lambda=lambdaROPE,alpha=alpha,target=I,cv=F)
ROPEKer = ROPEKernelEED(Y,lambda=lambdaROPE,alpha=alpha,target=I,cv=F)

AROPEPoisson = ROPEPoisson$IndROPEPoisson
AROPEKernel = ROPEKer$IndROPEKernel

#####################################################################################

# Method similar to ROPE, "RagsRidges" which utilizes partial correlation and Strimmer method:

OPT = optPenalty.LOOCV(Yvalid, lambdaMin = .5, lambdaMax = 30, step = 40,verbose = F,graph=F,target = I)

## Determine support regularized (standardized) precision under optimal penalty

ThetaRagsRidge = ridgeP(cov(Y),lambda=OPT$optLambda,target=I)

SparseThetaRagsRidge = sparsify(ThetaRagsRidge, threshold = "localFDR",verbose = F)

ARagsRidge = ifelse(SparseThetaRagsRidge$sparsePrecision != 0, 1, 0); diag(ARagsRidge) = 0

###################################################

PoissonG = graph.adjacency(AROPEPoisson,mode="undirected",diag=F)
ClustersPoisson = walktrap.community(PoissonG)

KernelG = graph.adjacency(AROPEKernel,mode="undirected",diag=F)
ClustersKernel = walktrap.community(KernelG)

RagsRidgeG = graph.adjacency(ARagsRidge,mode="undirected",diag=F)
ClustersRagsRidge = walktrap.community(RagsRidgeG)

#Plot graphs and clusters:
  
par(mfrow=c(2,2))

plot(ClustersTrue, TrueG, vertex.label=NA, vertex.size=5,main="Ground truth graph")
plot(ClustersPoisson, PoissonG, vertex.label=NA, vertex.size=5,main="ROPE (Poisson)")
plot(ClustersKernel, KernelG, vertex.label=NA, vertex.size=5,main="ROPE (Kernel)")
plot(ClustersRagsRidge, RagsRidgeG, vertex.label=NA, vertex.size=5,main="Rags2Ridges")
```

![graph1](https://user-images.githubusercontent.com/40263834/54021553-12a33b80-4199-11e9-98a5-1b0685fa46fa.png)


```r
# Larger sample size:

n = 200

X = mvrnorm(2*n,rep(0,p),Sigma)
X = huge.npn(X)

Yvalid = X[1:n,]
Y = X[(n+1):(2*n),]

lambdaROPE = ROPEPoissonEED(Yvalid,lambda=lambda,alpha=alpha,target=I)$lambda

ROPEPoisson = ROPEPoissonEED(Y,lambda=lambdaROPE,alpha=alpha,target=I,cv=F)
ROPEKer = ROPEKernelEED(Y,lambda=lambdaROPE,alpha=alpha,target=I,cv=F)

AROPEPoisson = ROPEPoisson$IndROPEPoisson
AROPEKernel = ROPEKer$IndROPEKernel

#####################################################################################

OPT = optPenalty.LOOCV(Yvalid, lambdaMin = .5, lambdaMax = 30, step = 40,verbose = F,graph=F,target = I)

ThetaRagsRidge = ridgeP(cov(Y),lambda=OPT$optLambda,target=I)

SparseThetaRagsRidge = sparsify(ThetaRagsRidge, threshold = "localFDR",verbose = F)

ARagsRidge = ifelse(SparseThetaRagsRidge$sparsePrecision != 0, 1, 0); diag(ARagsRidge) = 0

###################################################

PoissonG = graph.adjacency(AROPEPoisson,mode="undirected",diag=F)
ClustersPoisson = walktrap.community(PoissonG)

KernelG = graph.adjacency(AROPEKernel,mode="undirected",diag=F)
ClustersKernel = walktrap.community(KernelG)

RagsRidgeG = graph.adjacency(ARagsRidge,mode="undirected",diag=F)
ClustersRagsRidge = walktrap.community(RagsRidgeG)

par(mfrow=c(2,2))

plot(ClustersTrue, TrueG, vertex.label=NA, vertex.size=5,main="Ground truth graph")
plot(ClustersPoisson, PoissonG, vertex.label=NA, vertex.size=5,main="ROPE (Poisson)")
plot(ClustersKernel, KernelG, vertex.label=NA, vertex.size=5,main="ROPE (Kernel)")
plot(ClustersRagsRidge, RagsRidgeG, vertex.label=NA, vertex.size=5,main="Rags2Ridges")
```

![graph2](https://user-images.githubusercontent.com/40263834/54021962-079cdb00-419a-11e9-9a40-1d248c3fd6d2.png)

```r
# RagsRidge without FDR control:

# Calculate p-values from the RagsRdige Regularized partial correlations

R = pcor(ThetaRagsRidge) # Compute the partial correlation matrix

R = R[upper.tri(R)]

EmpiricDistribution = fdrtool(R,statistic = "correlation",verbose = F,plot=F) 

pvalues = EmpiricDistribution$pval # p-values determined from the partical correlation coefficients

ARagsRidgeNoFDR = matrix(0,p,p)

pvalues[pvalues > alpha] = 0 # Not significant p-values
pvalues[pvalues != 0] = 1

ARagsRidgeNoFDR[upper.tri(ARagsRidgeNoFDR)] = pvalues

ARagsRidgeNoFDR = ARagsRidgeNoFDR + t(ARagsRidgeNoFDR)

RagsRidgeGNoFDR = graph.adjacency(ARagsRidgeNoFDR,mode="undirected",diag=F)

ClustersRagsRidgeNoFDR = walktrap.community(RagsRidgeGNoFDR)

par(mfrow=c(1,2))

plot(ClustersRagsRidge, RagsRidgeG, vertex.label=NA, vertex.size=5,main="Rags2Ridges")
plot(ClustersRagsRidgeNoFDR, RagsRidgeGNoFDR, vertex.label=NA, vertex.size=5,main="Rags2Ridges (No FDR)")
```
![graph3](https://user-images.githubusercontent.com/40263834/54021977-12f00680-419a-11e9-8a77-ba24f88c3e6f.png)
