
# 7.5.2018

library(igraph)
library(spatstat)

SAureusNet = read.table("DREAM5_NetworkInference_GoldStandard_Network2.txt",sep="\t",header = F)

SAureusNet[1:10,]

SAureusNet = SAureusNet[,1:2]

colnames(SAureusNet) = c("col1","col2")

SAureusNet[1:10,]

IgraphSAureus = graph.data.frame(SAureusNet, directed = F)

Atrue = as_adjacency_matrix(IgraphSAureus,edges = F,type="both")

Atrue = as.matrix(Atrue)

isSymmetric(Atrue)

Atrue[1:5,1:5]

Atrue = ifelse(Atrue != 0, 1, 0)

dim(Atrue)

summary(colSums(Atrue))

Atrue[1:5,1:5]

f = function(m) t(m)[nrow(m):1,] # For plotting the images
plot(im(f(Atrue)), main="Real Adjacency",col = gray(seq(1,0,length.out=10)))

CW = cluster_walktrap(IgraphSAureus)

coords = layout_with_fr(IgraphSAureus, grid="nogrid")

length(table(membership(CW)))

NetColors = rainbow(length(table(membership(CW))))

GroupCol = NetColors[membership(CW)]


pdf("InSAureusNet.pdf", onefile=T, paper='A4r')

par(mfrow=c(1,1))

plot(CW, IgraphSAureus, mark.border=NA, mark.col = NA,layout=coords, vertex.label=NA, edge.arrow.size=.2,
     edge.color = "black", col=GroupCol, vertex.size = 2,main="In SAureus Network")

dev.off()

write.table(Atrue,"SAureusAdjacency.txt")
