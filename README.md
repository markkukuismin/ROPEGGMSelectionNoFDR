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
install.packages("FastGGM")

# Example

Explore network hubs and clusters from penalized precision matrix estimate:

```r
library(huge)
library(igraph)
library(rags2ridges)

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
```
Simulate data:

```r
v = 0.8
u = 0.1
L = huge.generator(n=n,d=p,g=4,graph="hub",v=v,u=u)  # Generate data with hub structures

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

ThetaRagsRidge = sparsify(ThetaRagsRidge, threshold = "localFDR",verbose = F)

ARagsRidge = ifelse(ThetaRagsRidge$sparsePrecision != 0, 1, 0); diag(ARagsRidge) = 0

###################################################

PoissonG = graph.adjacency(AROPEPoisson,mode="undirected",diag=F)
KernelG = graph.adjacency(AROPEKernel,mode="undirected",diag=F)
RagsRidgeG = graph.adjacency(ARagsRidge,mode="undirected",diag=F)
```
Plot graphs:

```r
par(mfrow=c(2,2))

plot(TrueG, vertex.label=NA, vertex.size=5,main="Ground truth graph")
plot(PoissonG, vertex.label=NA, vertex.size=5,main="ROPE (Poisson)")
plot(KernelG, vertex.label=NA, vertex.size=5,main="ROPE (Kernel)")
plot(RagsRidgeG, vertex.label=NA, vertex.size=5,main="Rags2Ridges")
```
![graphfigures](https://user-images.githubusercontent.com/40263834/51476596-a5754b80-1d8e-11e9-9fe7-7f6e9571608f.png)

Change n = 200 (400 samples) and five (5) hubs:

![graphfigures2](https://user-images.githubusercontent.com/40263834/51476807-495ef700-1d8f-11e9-9842-b7f8e2ed9c10.png)
