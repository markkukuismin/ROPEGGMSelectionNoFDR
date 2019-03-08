
# 30.5.2018


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

#S = cor(Y)
#>summary(S[upper.tri(S)])
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.94482 -0.10929  0.08497  0.11492  0.35207  0.99001 
#>summary(abs(S[upper.tri(S)]))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.09633 0.23197 0.27933 0.42183 0.99001 

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

# How FastGGM finds:

# 1) Transciption factors. Are they identified as hubs?

# 2) Pathogen cluster

#####################################################

FGGM = FastGGM_Parallel(Y)

AFGGM = ifelse(FGGM$CI_low_parCor*FGGM$CI_high_parCor < 0, 0, 1)
diag(AFGGM) = 0

colnames(AFGGM) = GeneNames$Names

FGGMGraph = graph.adjacency(AFGGM, mode="undirected")

DegreeFGGM = colSums(AFGGM)

# 1) Are transcription factors identified as hubs?

# Examine, say, 100 hub nodes:

NmbSelected = 100

HubsFGGM = order(DegreeFGGM,decreasing = T)[1:NmbSelected]

HubsFGGM = GeneNames$Names[HubsFGGM]

length(intersect(HubsFGGM,TranscriptionFactors))

intersect(HubsFGGM,TranscriptionFactors)

as.character(HubsFGGM[1:10])

nmb = p*(p-1)/2

100*sum(AFGGM[upper.tri(AFGGM)])/nmb

###############################################

ClustersFGGM = walktrap.community(FGGMGraph)

length(table(membership(ClustersFGGM)))


#########

coordsFGGM = layout_with_fr(FGGMGraph, grid="nogrid")

ESizes = rep(2,p)

V(FGGMGraph)$size=ESizes

NetColors = rep("white",p)

PathogenClusterInd = rep(NA,27)

for(i in 1:27){
  
  PathogenClusterInd[i] = which(PathogenCluster[i] == GeneNames$Names) 
  
}

TranscriptionInd = 1:99

NetColors[PathogenClusterInd] = "lightgreen"
NetColors[TranscriptionInd] = "red"


plot(ClustersFGGM, FGGMGraph, mark.border=NA, mark.col = NA,layout=coordsFGGM, vertex.label=NA,
     edge.arrow.size=.2, edge.color = "black",main="FastGGM graph",col=NetColors)


##################################################################

# 2) How is the cluster enrichted for pathogenic genes identified?

table(ClustersFGGM$membership[PathogenClusterInd])


AFGGMPathogen = AFGGM[PathogenClusterInd,PathogenClusterInd]

FGGMGraphPathogen = graph.adjacency(AFGGMPathogen,mode="undirected")

ClustersFGGMPathogen = walktrap.community(FGGMGraphPathogen)

coordsFGGM = layout_with_fr(FGGMGraphPathogen, grid="nogrid",niter=10000)

table(ClustersFGGM$membership[PathogenClusterInd])
PathogenColsFGGM = rainbow(
  max(ClustersFGGM$membership[PathogenClusterInd]))[ClustersFGGM$membership[PathogenClusterInd]]

FGGMPathodenDergree = colSums(AFGGMPathogen)
FGGMPathodenNodeSize = 200*(FGGMPathodenDergree/sum(FGGMPathodenDergree))

V(FGGMGraphPathogen)$label.cex = 0.7
V(FGGMGraphPathogen)$size = FGGMPathodenNodeSize

PathogenShapes = rep("circle",27)
PathogenShapes[c(1,4)] = "square"

pdf("FastGGMPathogenGraph.pdf", onefile=T, paper='A4r')

plot(ClustersFGGMPathogen, FGGMGraphPathogen, mark.border=NA, mark.col = NA,layout=coordsFGGM, 
     vertex.label=PathogenCluster, edge.arrow.size=.2, edge.color = "gray",
     main="FastGGM Pathogen network",col=PathogenColsFGGM, vertex.shape=PathogenShapes)

dev.off()
