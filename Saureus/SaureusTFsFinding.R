
# 22.05.2018

library(c3net)
library(rags2ridges)
library(camel)
library(BDgraph)
library(GGMridge)
library(huge)

source("GGMRidge.txt")
source("SpSeFallPreMCCVector.txt")
source("ROPE.txt")
source("ROPEEEDEmpiricPoisson.r")
source("ROPEEEDKernelMethod.r")
source("Scale.r")
source("crossvalidationROPE.txt")

Data = read.table("net2_expression_data.tsv",header=T)

GeneNames = read.table("net2_gene_ids.tsv",skip=1,header=F)

TranscriptionFactors = read.table("net2_transcription_factors.tsv",header = F)

colnames(GeneNames) = c("ID","Names")
colnames(Data) = GeneNames[,2]

Data[1:3,1:5]
GeneNames[1:5,]

# Transcription factors are in order:

TranscriptionFactors =  GeneNames$Names[1:99]

dim(Data)

# From DREAM5 web page:

# "The expression data has been uniformly normalized for each compendium so 
#  that values are comparable across experiments."

Y = as.matrix(Data)

# However, many methods in the Marbach et al (2012) standardized gene profiles
# to zero mean and unit variance.

#Y = huge.npn(Data)
Y = scale(Y)

p = ncol(Y)
n = nrow(Y)

hist(Y[,sample(1:p,1)],probability = T)

# How methods identify hubs and how they match transcription factors?

# According to Marbach et al (2012) these are pathogens form a module of their own:

PathogenCluster = c("SAV2357","SAV0814","SAV2503","SAV2553","SAV0825","SAV1049","SAV0717",
                    "SAV0715","SAV0718","SAV0708","SAV0716","SAV0849","SAV0435","SAV1431",
                    "SAV0714","SAV0423","SAV0429","SAV2221","SAV1105","SAV0220","SAV1159",
                    "SAV0229","SAV2001","SAV2222","SAV0424" ,"SAV0426","SAV0749")

# These are Transcription factors in the pathogen cluster:

PathogenClusterTranscriptions =  intersect(PathogenCluster,TranscriptionFactors)

#####################################################

# How different methods find:

# 1) Transciption factors. Are they identified as hubs?

# 2) Pathogen cluster

#####################################################

alpha = 0.1

I = diag(1,p)

# ROPE

#lambda = seq(1,100,length.out = 20)
#ROPEPoisson = ROPEPoissonEED(Y, lambda=lambda, alpha=alpha,cv=T, target=I)
#ROPEPoisson$lambda
#lambdaROPE = ROPEPoisson$lambda
#plot(ROPEPoisson$cv.values)

lambda = 5

ROPEPoisson = ROPEPoissonEED(Y, lambda=lambda, alpha=alpha, cv=F, target=I)
ROPEKernel = ROPEKernelEED(Y, lambda=lambda, alpha=alpha, cv=F, target=I)

AROPEPoisson = ROPEPoisson$IndROPE
AROPEKernel = ROPEKernel$IndROPE

diag(AROPEPoisson) = 0
diag(AROPEKernel) = 0

colnames(AROPEPoisson) = GeneNames$Names
colnames(AROPEKernel) = GeneNames$Names

ROPEPoissonGraph = graph.adjacency(AROPEPoisson, mode="undirected")
ROPEKernelGraph = graph.adjacency(AROPEKernel, mode="undirected")

DegreeROPEPoisson = colSums(AROPEPoisson)
DegreeROPEKernel = colSums(AROPEKernel)

# 1) Are transcription factors identified as hubs?

# Examine, say, 100 hub nodes:

NmbSelected = 100

HubsPoisson = order(DegreeROPEPoisson,decreasing = T)[1:NmbSelected]
HubsKernel = order(DegreeROPEKernel,decreasing = T)[1:NmbSelected]

HubsPoisson = GeneNames$Names[HubsPoisson]
HubsKernel = GeneNames$Names[HubsKernel]

length(intersect(HubsPoisson,TranscriptionFactors))
length(intersect(HubsKernel,TranscriptionFactors))

intersect(HubsPoisson,TranscriptionFactors)

intersect(HubsKernel,TranscriptionFactors)

as.character(HubsPoisson[1:10])

# 2) How is the cluster enrichted for pathogenic genes identified?

ClustersPoisson = walktrap.community(ROPEPoissonGraph)
ClustersKernel = walktrap.community(ROPEKernelGraph)

length(table(membership(ClustersPoisson)))

length(table(membership(ClustersKernel)))

#########

coords = layout_with_fr(ROPEPoissonGraph, grid="nogrid")

ESizes = rep(2,p)

V(ROPEPoissonGraph)$size=ESizes

NetColors = rep("white",p)

PathogenClusterInd = rep(NA,27)

for(i in 1:27){
  
  PathogenClusterInd[i] = which(PathogenCluster[i] == GeneNames$Names) 
  
}

TranscriptionInd = 1:99

NetColors[PathogenClusterInd] = "lightgreen"
NetColors[TranscriptionInd] = "red"

plot(ClustersPoisson, ROPEPoissonGraph, mark.border=NA, mark.col = NA,layout=coords, vertex.label=NA,
     edge.arrow.size=.2, edge.color = "black",main="ROPE Poisson graph",col=NetColors)

####

coords = layout_with_fr(ROPEKernelGraph, grid="nogrid")

V(ROPEKernelGraph)$size=ESizes

plot(ClustersKernel, ROPEKernelGraph, mark.border=NA, mark.col = NA,layout=coords, vertex.label=NA,
     edge.arrow.size=.2, edge.color = "black",main="ROPE Kernel graph",col=NetColors)

#####################################################
#####################################################

# RagsRidge


#lambdaRagsRidge  = optPenalty.aLOOCV(Y, lambdaMin = 0.01, lambdaMax = 5, step = 20, 
#                                     verbose = T, type = "Alt", target = I, graph=F, output="all")

#lambdaRagsRidge = lambdaRagsRidge$optLambda # Penalty chosen with leave-one-out crossvalidation

#> lambdaRagsRidge
#[1] 0.01

lambdaRagsRidge = 5

ThetaRR = ridgeP(cov(Y), lambda=lambdaRagsRidge, type = "Alt", target = I)

ThetaRagsRidge = sparsify(ThetaRR, threshold = "localFDR",verbose = F)

ThetaRagsRidge = ThetaRagsRidge$sparsePrecision

ARRidge = ifelse(ThetaRagsRidge != 0, 1, 0); diag(ARRidge) = 0

RagsRidgeGraph = graph.adjacency(ARRidge, mode="undirected")

DegreeRagsRidge = colSums(ARRidge)

# 1) Are transcription factors identified as hubs?

# Examine, say, 100 hub nodes:

HubsRagsRidge = order(DegreeRagsRidge,decreasing = T)[1:NmbSelected]

HubsRagsRidge = GeneNames$Names[HubsRagsRidge]

length(intersect(HubsRagsRidge,TranscriptionFactors))

intersect(HubsRagsRidge,TranscriptionFactors)


# 2) How is the cluster enrichted for pathogenic genes identified?

ClustersRagsRidge = walktrap.community(RagsRidgeGraph)

length(table(membership(ClustersRagsRidge)))


#########

coords = layout_with_fr(RagsRidgeGraph, grid="nogrid")

V(RagsRidgeGraph)$size=ESizes

plot(ClustersRagsRidge, RagsRidgeGraph, mark.border=NA, mark.col = NA,layout=coords, vertex.label=NA,
     edge.arrow.size=.2, edge.color = "black",main="RagsRidge graph",col=NetColors)

#########################################################
#########################################################

# GGMRidge

GGMRidgeEst = GGMRidge(Y,cv=T)

RRidge = GGMRidgeEst$GGMRidgePartialCor

RRidge[is.na(RRidge)] = 0

AGGMRidge = ifelse(RRidge != 0, 1, 0); diag(AGGMRidge) = 0

GGMRidgeGraph = graph.adjacency(AGGMRidge, mode="undirected")

DegreeGGMRidge = colSums(AGGMRidge)

# 1) Are transcription factors identified as hubs?

# Examine, say, 100 hub nodes:

HubsGGMRidge = order(DegreeGGMRidge,decreasing = T)[1:NmbSelected]

HubsGGMRidge = GeneNames$Names[HubsGGMRidge]

length(intersect(HubsGGMRidge,TranscriptionFactors))

intersect(HubsGGMRidge,TranscriptionFactors)


# 2) How is the cluster enrichted for pathogenic genes identified?

ClustersGGMRidge = walktrap.community(GGMRidgeGraph)

length(table(membership(ClustersGGMRidge)))


#########

coords = layout_with_fr(GGMRidgeGraph, grid="nogrid")

V(GGMRidgeGraph)$size=ESizes

plot(ClustersGGMRidge, GGMRidgeGraph, mark.border=NA, mark.col = NA,layout=coords, vertex.label=NA,
     edge.arrow.size=.2, edge.color = "black",main="GGMRidge graph",col=NetColors)


#####################################################

# TIGER


OutTIGER = camel.tiger(Y,method="slasso",verbose=F)
TIGERResults = camel.tiger.select(OutTIGER,criterion="stars",stars.thresh = 0.05,verbose=T)

TIGERResults = camel.tiger(Y,method="slasso",verbose=F,lambda=TIGERResults$opt.lambda)

ThetaTIGER = TIGERResults$icov
ThetaTIGER = as.matrix(ThetaTIGER[[1]])

ATIGER = ifelse(ThetaTIGER != 0, 1, 0)
diag(ATIGER) = 0
TIGERGraph = graph.adjacency(ATIGER, mode="undirected")

DegreeTIGER = colSums(ATIGER)

# 1) Are transcription factors identified as hubs?

# Examine, say, 100 hub nodes:

HubsTIGER = order(DegreeTIGER,decreasing = T)[1:NmbSelected]

HubsTIGER = GeneNames$Names[HubsTIGER]

length(intersect(HubsTIGER,TranscriptionFactors))

intersect(HubsTIGER,TranscriptionFactors)


# 2) How is the cluster enrichted for pathogenic genes identified?

ClustersTIGER = walktrap.community(TIGERGraph)

length(table(membership(ClustersTIGER)))


#########

coords = layout_with_fr(TIGERGraph, grid="nogrid")

V(TIGERGraph)$size=ESizes

plot(ClustersTIGER, TIGERGraph, mark.border=NA, mark.col = NA,layout=coords, vertex.label=NA,
     edge.arrow.size=.2, edge.color = "black",main="TIGER graph",col=NetColors)

#######################################
#######################################

# A better look at the pathogen genes:

# 2) How is the cluster enrichted for pathogenic genes identified?

table(ClustersPoisson$membership[PathogenClusterInd])
table(ClustersKernel$membership[PathogenClusterInd])

table(ClustersRagsRidge$membership[PathogenClusterInd])
table(ClustersGGMRidge$membership[PathogenClusterInd])

table(ClustersTIGER$membership[PathogenClusterInd])


########################################

APoissonPathogen = AROPEPoisson[PathogenClusterInd,PathogenClusterInd]

ROPEPoissonGraphPathogen = graph.adjacency(APoissonPathogen,mode="undirected")

ClustersPoissonPathogen = walktrap.community(ROPEPoissonGraphPathogen)

coordsRP = layout_with_fr(ROPEPoissonGraphPathogen, grid="nogrid",niter=10000)

table(ClustersPoisson$membership[PathogenClusterInd])
PathogenColsRP = rainbow(
  max(ClustersPoisson$membership[PathogenClusterInd]))[ClustersPoisson$membership[PathogenClusterInd]]

PathogenShapes = rep("circle",27)
PathogenShapes[c(1,4)] = "square"

RPPathodenDergree = colSums(APoissonPathogen)
RPPathodenNodeSize = 200*(RPPathodenDergree/sum(RPPathodenDergree))

V(ROPEPoissonGraphPathogen)$label.cex = 0.7
V(ROPEPoissonGraphPathogen)$size = RPPathodenNodeSize

pdf("ROPEPoissonPathogenGraph.pdf", onefile=T, paper='A4r')

plot(ClustersPoissonPathogen, ROPEPoissonGraphPathogen, mark.border=NA, mark.col = NA,layout=coordsRP, 
     vertex.label=PathogenCluster, edge.arrow.size=.2, edge.color = "gray",
     main="ROPE Poisson Pathogen network",col=PathogenColsRP, vertex.shape=PathogenShapes)

dev.off()

########################################

AKernelPathogen = AROPEKernel[PathogenClusterInd,PathogenClusterInd]

ROPEKernelGraphPathogen = graph.adjacency(AKernelPathogen,mode="undirected")

ClustersKernelPathogen = walktrap.community(ROPEKernelGraphPathogen)

coordsRK = layout_with_fr(ROPEKernelGraphPathogen, grid="nogrid",niter=10000)

table(ClustersKernel$membership[PathogenClusterInd])
PathogenColsRK = rainbow(
  max(ClustersKernel$membership[PathogenClusterInd]))[ClustersKernel$membership[PathogenClusterInd]]

RKPathodenDergree = colSums(AKernelPathogen)
RKPathodenNodeSize = 200*(RKPathodenDergree/sum(RKPathodenDergree))

V(ROPEKernelGraphPathogen)$label.cex = 0.7
V(ROPEKernelGraphPathogen)$size = RKPathodenNodeSize

pdf("ROPEKernelPathogenGraph.pdf", onefile=T, paper='A4r')

plot(ClustersKernelPathogen, ROPEKernelGraphPathogen, mark.border=NA, mark.col = NA,layout=coordsRK, 
     vertex.label=PathogenCluster, edge.arrow.size=.2, edge.color = "gray",
     main="ROPE Kernel Pathogen network",col=PathogenColsRK, vertex.shape=PathogenShapes)

dev.off()

########################################

ARagsRidgePathogen = ARRidge[PathogenClusterInd,PathogenClusterInd]

RagsRidgeGraphPathogen = graph.adjacency(ARagsRidgePathogen,mode="undirected")

ClustersRagsRidgePathogen = walktrap.community(RagsRidgeGraphPathogen)

coordsRR = layout_with_fr(RagsRidgeGraphPathogen, grid="nogrid",niter=10000)

table(ClustersRagsRidge$membership[PathogenClusterInd])
PathogenColsRR = rainbow(
  max(ClustersRagsRidge$membership[PathogenClusterInd]))[ClustersRagsRidge$membership[PathogenClusterInd]]

plot(ClustersRagsRidgePathogen, RagsRidgeGraphPathogen, mark.border=NA, mark.col = NA,layout=coordsRR, 
     vertex.label=PathogenCluster, edge.arrow.size=.2, edge.color = "gray",
     main="RagsRidge Pathogen network",col=PathogenColsRR, vertex.shape=PathogenShapes)

########################################

AGGMRidgePathogen = AGGMRidge[PathogenClusterInd,PathogenClusterInd]

GGMRidgeGraphPathogen = graph.adjacency(AGGMRidgePathogen,mode="undirected")

ClustersGGMRidgePathogen = walktrap.community(GGMRidgeGraphPathogen)

coordsGGM = layout_with_fr(GGMRidgeGraphPathogen, grid="nogrid",niter=10000)

table(ClustersGGMRidge$membership[PathogenClusterInd])
PathogenColsGGM = rainbow(
  max(ClustersGGMRidge$membership[PathogenClusterInd]))[ClustersGGMRidge$membership[PathogenClusterInd]]

plot(ClustersGGMRidgePathogen, GGMRidgeGraphPathogen, mark.border=NA, mark.col = NA,layout=coordsGGM, 
     vertex.label=PathogenCluster, edge.arrow.size=.2, edge.color = "gray",
     main="GGMRidge Pathogen network",col=PathogenColsGGM, vertex.shape=PathogenShapes)

########################################

ATIGERPathogen = ATIGER[PathogenClusterInd,PathogenClusterInd]

TIGERGraphPathogen = graph.adjacency(ATIGERPathogen,mode="undirected")

ClustersTIGERPathogen = walktrap.community(TIGERGraphPathogen)

coordsT = layout_with_fr(TIGERGraphPathogen, grid="nogrid",niter=10000)

table(ClustersTIGER$membership[PathogenClusterInd])
PathogenColsT = rainbow(
  max(ClustersTIGER$membership[PathogenClusterInd]))[ClustersTIGER$membership[PathogenClusterInd]]

TIGERPathodenDergree = colSums(ATIGERPathogen)
TIGERPathodenNodeSize = 200*(TIGERPathodenDergree/sum(TIGERPathodenDergree))

V(TIGERGraphPathogen)$label.cex = 0.7
V(TIGERGraphPathogen)$size = TIGERPathodenNodeSize

pdf("TIGERPathogenGraph.pdf", onefile=T, paper='A4r')

plot(ClustersTIGERPathogen, TIGERGraphPathogen, mark.border=NA, mark.col = NA,layout=coordsT, 
     vertex.label=PathogenCluster, edge.arrow.size=.2, edge.color = "gray",
     main="TIGER Pathogen network",col=PathogenColsT, vertex.shape=PathogenShapes)

dev.off()

########################################

nmb = p*(p-1)/2

100*sum(AROPEPoisson[upper.tri(AROPEPoisson)])/nmb
100*sum(AROPEKernel[upper.tri(AROPEKernel)])/nmb

100*sum(ARRidge[upper.tri(ARRidge)])/nmb
100*sum(AGGMRidge[upper.tri(AGGMRidge)])/nmb

100*sum(ATIGER[upper.tri(ATIGER)])/nmb

########################################

# Let's look at the hubs once moore:

as.character(HubsPoisson[1:10])
as.character(HubsKernel[1:10])

as.character(HubsRagsRidge[1:10])
as.character(HubsGGMRidge[1:10])

as.character(HubsTIGER[1:10])

intersect(HubsPoisson[1:10],HubsKernel[1:10])

intersect(HubsPoisson[1:10],HubsTIGER[1:10])

which(GeneNames$Names == "SAV2011")


SAV2011NeighboursPoisson = which(AROPEPoisson[,2602] != 0)
SAV2011DataP = data.frame(ID=GeneNames$ID[SAV2011NeighboursPoisson],
                          Name=GeneNames$Names[SAV2011NeighboursPoisson])
write.table(SAV2011DataP,"SAV2011NeighboursPoisson.txt",row.names = F,quote=F)

SAV2011NeighboursKernel = which(AROPEKernel[,2602] != 0)
SAV2011DataK = data.frame(ID=GeneNames$ID[SAV2011NeighboursKernel],Name=GeneNames$Names[SAV2011NeighboursKernel])
write.table(SAV2011DataK,"SAV2011NeighboursKernel.txt",row.names = F,quote=F)

membership(ClustersPoisson)[2602]

SAV2011Cluster = which(membership(ClustersPoisson) == 1)

SAV2011Adjacency = AROPEPoisson[SAV2011Cluster,SAV2011Cluster]

SAV2011Graph = graph.adjacency(SAV2011Adjacency,mode="undirected")

SAV2011WTC = walktrap.community(SAV2011Graph)

coords = layout_with_fr(SAV2011Graph, grid="nogrid",niter=10000)

# This is also the largest cluster:

plot(SAV2011WTC, SAV2011Graph, mark.border=NA, mark.col = NA,layout=coords, 
     vertex.label=NA, edge.arrow.size=.2, edge.color = "black",
     main="Cluster of the highest hub node",col="white")
