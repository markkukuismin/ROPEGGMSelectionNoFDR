
# 04.03.2019

library(ggplot2)
library(reshape)
library(grid)
library(gridExtra)
library(scales)
library(gtable)

Model = "Hub"
Dimension = 500
Samplesize = c(floor(0.1*Dimension),floor(0.25*Dimension),floor(0.5*Dimension),floor(1*Dimension))

Data1 = read.csv(paste(Model,"NmbResultsForDimension",Dimension,"n=",Samplesize[1],".csv",sep=""),header=T)
Data2 = read.csv(paste(Model,"NmbResultsForDimension",Dimension,"n=",Samplesize[2],".csv",sep=""),header=T)
Data3 = read.csv(paste(Model,"NmbResultsForDimension",Dimension,"n=",Samplesize[3],".csv",sep=""),header=T)
Data4 = read.csv(paste(Model,"NmbResultsForDimension",Dimension,"n=",Samplesize[4],".csv",sep=""),header=T)

TrueNmb = Data1$TrueNmbOfCluster[1]

Names = colnames(Data1)

Names[Names=="TIGERStARS"] = "TIGER"
Names[Names=="TrueNmbOfClusters"] = "p"

colnames(Data1) = Names
colnames(Data2) = Names
colnames(Data3) = Names
colnames(Data4) = Names

rm(Names)

MeltData1 = melt(Data1, id=c('p'))
MeltData2 = melt(Data2, id=c('p'))
MeltData3 = melt(Data3, id=c('p'))
MeltData4 = melt(Data4, id=c('p'))

colnames(MeltData1) = colnames(MeltData2) = colnames(MeltData3) = colnames(MeltData4) = c("p","Method","value")

p1 = ggplot(MeltData1) +
  geom_boxplot(aes(x=Method, y=value, color=Method)) +
  facet_grid(Dimension) + 
  theme(axis.text=element_text(size=8),axis.title=element_text(size=10,face="bold")) + 
  ggtitle(paste(Model," graph, n=",Samplesize[1],sep="")) +
  labs(y="Number of identified clusters",x=NULL) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_text(size=8),
        legend.text=element_text(size=6)) + 
  geom_hline(yintercept = TrueNmb, linetype=2)

p2 = ggplot(MeltData2) +
  geom_boxplot(aes(x=Method, y=value, color=Method)) +
  facet_grid(Dimension) + 
  theme(axis.text=element_text(size=8),axis.title=element_text(size=10,face="bold")) + 
  ggtitle(paste(Model," graph, n=",Samplesize[2],sep="")) +
  labs(y="Number of identified clusters",x=NULL) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_text(size=8),
        legend.text=element_text(size=6)) + 
  geom_hline(yintercept = TrueNmb, linetype=2)

p3 = ggplot(MeltData3) +
  geom_boxplot(aes(x=Method, y=value, color=Method)) +
  facet_grid(Dimension) + 
  theme(axis.text=element_text(size=8),axis.title=element_text(size=10,face="bold")) + 
  ggtitle(paste(Model," graph, n=",Samplesize[3],sep="")) +
  labs(y="Number of identified clusters",x=NULL) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_text(size=8),
        legend.text=element_text(size=6)) + 
  geom_hline(yintercept = TrueNmb, linetype=2)

p4 = ggplot(MeltData4) +
  geom_boxplot(aes(x=Method, y=value, color=Method)) +
  facet_grid(Dimension) + 
  theme(axis.text=element_text(size=8),axis.title=element_text(size=10,face="bold")) + 
  ggtitle(paste(Model," graph, n=",Samplesize[4],sep="")) +
  labs(y="Number of identified clusters",x=NULL) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_text(size=8),
        legend.text=element_text(size=6)) + 
  geom_hline(yintercept = TrueNmb, linetype=2)

############################################################

grid.arrange(p1, p2, p3, p4, ncol=2)

############################################################

scale=0.8

pdf(paste(Dimension,Model,"FoundClusters.pdf",sep=""), width = 10*scale, height = 7*scale)
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off() 

############################################################

if(Model == "Cluster"){
  
  Data1 = read.csv(paste(Model,"WholeGraphDiagnosticForDimension",Dimension,"n=",Samplesize[1],".csv",sep=""),header=T)
  Data2 = read.csv(paste(Model,"WholeGraphDiagnosticForDimension",Dimension,"n=",Samplesize[2],".csv",sep=""),header=T)
  Data3 = read.csv(paste(Model,"WholeGraphDiagnosticForDimension",Dimension,"n=",Samplesize[3],".csv",sep=""),header=T)
  Data4 = read.csv(paste(Model,"WholeGraphDiagnosticForDimension",Dimension,"n=",Samplesize[4],".csv",sep=""),header=T)
  
  Data1 = Data1[,c(2,4,5,6)]
  Data2 = Data2[,c(2,4,5,6)]
  Data3 = Data3[,c(2,4,5,6)]
  Data4 = Data4[,c(2,4,5,6)]
  
  Data1[is.na(Data1)] = 0
  Data2[is.na(Data2)] = 0
  Data3[is.na(Data3)] = 0
  Data4[is.na(Data4)] = 0
  
  MeltData1 = melt(Data1, id=c('Method'))
  MeltData2 = melt(Data2, id=c('Method'))
  MeltData3 = melt(Data3, id=c('Method'))
  MeltData4 = melt(Data4, id=c('Method'))
  
  
  p1 = ggplot(MeltData1) +
    geom_boxplot(aes(x=Method, y=value, color=Method))+
    facet_grid(.~variable) + 
    theme(axis.text=element_text(size=8),axis.title=element_text(size=10,face="bold")) + 
    ggtitle(paste(Model," graph, p=",Dimension, ", n=",Samplesize[1],sep="")) +
    labs(x=NULL) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(color = "black", size=8, vjust=.8, hjust=0.8),
          axis.ticks.x=element_blank(),
          legend.title=element_text(size=5),
          legend.text=element_text(size=4),
          legend.key.size = unit(1,"line")) +
    scale_x_discrete(labels=NULL)
  
  p2 = ggplot(MeltData2) +
    geom_boxplot(aes(x=Method, y=value, color=Method))+
    facet_grid(.~variable) + 
    theme(axis.text=element_text(size=8),axis.title=element_text(size=10,face="bold")) + 
    ggtitle(paste(Model," graph, p=",Dimension, ", n=",Samplesize[2],sep="")) +
    labs(x=NULL) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(color = "black", size=8, vjust=.8, hjust=0.8),
          axis.ticks.x=element_blank(),
          legend.title=element_text(size=5),
          legend.text=element_text(size=4),
          legend.key.size = unit(1,"line")) +
    scale_x_discrete(labels=NULL)
  
  p3 = ggplot(MeltData3) +
    geom_boxplot(aes(x=Method, y=value, color=Method))+
    facet_grid(.~variable) + 
    theme(axis.text=element_text(size=8),axis.title=element_text(size=10,face="bold")) + 
    ggtitle(paste(Model," graph, p=",Dimension, ", n=",Samplesize[3],sep="")) +
    labs(x=NULL) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(color = "black", size=8, vjust=.8, hjust=0.8),
          axis.ticks.x=element_blank(),
          legend.title=element_text(size=5),
          legend.text=element_text(size=4),
          legend.key.size = unit(1,"line")) +
    scale_x_discrete(labels=NULL)
  
  p4 = ggplot(MeltData4) +
    geom_boxplot(aes(x=Method, y=value, color=Method))+
    facet_grid(.~variable) + 
    theme(axis.text=element_text(size=8),axis.title=element_text(size=10,face="bold")) + 
    ggtitle(paste(Model," graph, p=",Dimension, ", n=",Samplesize[4],sep="")) +
    labs(x=NULL) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(color = "black", size=8, vjust=.8, hjust=0.8),
          axis.ticks.x=element_blank(),
          legend.title=element_text(size=5),
          legend.text=element_text(size=4),
          legend.key.size = unit(1,"line")) +
    scale_x_discrete(labels=NULL)
  
  
  ############################################################
  
  grid.arrange(p1, p2, p3, p4, ncol=2)
  
  ############################################################
  
  scale=0.8
  
  pdf(paste(Dimension,Model,"WholeGraphDiagnostics.pdf",sep=""), width = 11.5*scale, height = 8*scale)
  grid.arrange(p1, p2, p3, p4, ncol=2)
  dev.off() 
  
  ############################################################
  
}

if(Model == "Hub"){
  
  Data1 = read.csv(paste(Model,"NmbOfCorretlyIdentHubsForDimension",Dimension,"n=",Samplesize[1],".csv",sep=""),header=T)
  Data2 = read.csv(paste(Model,"NmbOfCorretlyIdentHubsForDimension",Dimension,"n=",Samplesize[2],".csv",sep=""),header=T)
  Data3 = read.csv(paste(Model,"NmbOfCorretlyIdentHubsForDimension",Dimension,"n=",Samplesize[3],".csv",sep=""),header=T)
  Data4 = read.csv(paste(Model,"NmbOfCorretlyIdentHubsForDimension",Dimension,"n=",Samplesize[4],".csv",sep=""),header=T)
  
  MeltData1 = melt(Data1)
  MeltData2 = melt(Data2)
  MeltData3 = melt(Data3)
  MeltData4 = melt(Data4)
  
  colnames(MeltData1) = colnames(MeltData2) = colnames(MeltData3) = colnames(MeltData4) = c("Method","Nmb")
  
  p1 = ggplot(MeltData1) +
    geom_boxplot(aes(x=Method, y=Nmb, color=Method)) +
    facet_grid(Dimension) + 
    theme(axis.text=element_text(size=8),axis.title=element_text(size=10,face="bold")) + 
    ggtitle(paste(Model," graph, n=",Samplesize[1],sep="")) +
    labs(y="Number of identified hubs",x=NULL) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_blank(),
          axis.ticks.x=element_blank(),
          legend.title=element_text(size=8),
          legend.text=element_text(size=6)) + 
    geom_hline(yintercept = TrueNmb, linetype=2)
  
  p2 = ggplot(MeltData2) +
    geom_boxplot(aes(x=Method, y=Nmb, color=Method)) +
    facet_grid(Dimension) + 
    theme(axis.text=element_text(size=8),axis.title=element_text(size=10,face="bold")) + 
    ggtitle(paste(Model," graph, n=",Samplesize[2],sep="")) +
    labs(y="Number of identified hubs",x=NULL) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_blank(),
          axis.ticks.x=element_blank(),
          legend.title=element_text(size=8),
          legend.text=element_text(size=6)) + 
    geom_hline(yintercept = TrueNmb, linetype=2)
  
  p3 = ggplot(MeltData3) +
    geom_boxplot(aes(x=Method, y=Nmb, color=Method)) +
    facet_grid(Dimension) + 
    theme(axis.text=element_text(size=8),axis.title=element_text(size=10,face="bold")) + 
    ggtitle(paste(Model," graph, n=",Samplesize[3],sep="")) +
    labs(y="Number of identified hubs",x=NULL) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_blank(),
          axis.ticks.x=element_blank(),
          legend.title=element_text(size=8),
          legend.text=element_text(size=6)) + 
    geom_hline(yintercept = TrueNmb, linetype=2)
  
  p4 = ggplot(MeltData4) +
    geom_boxplot(aes(x=Method, y=Nmb, color=Method)) +
    facet_grid(Dimension) + 
    theme(axis.text=element_text(size=8),axis.title=element_text(size=10,face="bold")) + 
    ggtitle(paste(Model," graph, n=",Samplesize[4],sep="")) +
    labs(y="Number of identified hubs",x=NULL) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_blank(),
          axis.ticks.x=element_blank(),
          legend.title=element_text(size=8),
          legend.text=element_text(size=6)) + 
    geom_hline(yintercept = TrueNmb, linetype=2)
  
  ############################################################
  
  grid.arrange(p1, p2, p3, p4, ncol=2)
  
  ############################################################
  
  scale=0.8
  
  pdf(paste(Dimension,Model,"NumberOfIdentifiedHubs.pdf",sep=""), width = 10*scale, height = 7*scale)
  grid.arrange(p1, p2, p3, p4, ncol=2)
  dev.off() 
  
  ############################################################
  
  Neigh = c("First","Second","Third")
  
  for(i in 1:TrueNmb){
    
    Data1 = read.csv(paste(Model,Neigh[i],"NeighborhoodDiagnosticsForDimension",Dimension,"n=",Samplesize[1],
                           ".csv",sep=""),header=T,row.names=1)
    Data2 = read.csv(paste(Model,Neigh[i],"NeighborhoodDiagnosticsForDimension",Dimension,"n=",Samplesize[2],
                           ".csv",sep=""),header=T,row.names=1)
    Data3 = read.csv(paste(Model,Neigh[i],"NeighborhoodDiagnosticsForDimension",Dimension,"n=",Samplesize[3],
                           ".csv",sep=""),header=T,row.names=1)
    Data4 = read.csv(paste(Model,Neigh[i],"NeighborhoodDiagnosticsForDimension",Dimension,"n=",Samplesize[4],
                           ".csv",sep=""),header=T,row.names=1)
    
    Data1 = Data1[,c(2,4,5,6)]
    Data2 = Data2[,c(2,4,5,6)]
    Data3 = Data3[,c(2,4,5,6)]
    Data4 = Data4[,c(2,4,5,6)]
    
    Data1[is.na(Data1)] = 0
    Data2[is.na(Data2)] = 0
    Data3[is.na(Data3)] = 0
    Data4[is.na(Data4)] = 0
    
    MeltData1 = melt(Data1, id=c('Method'))
    MeltData2 = melt(Data2, id=c('Method'))
    MeltData3 = melt(Data3, id=c('Method'))
    MeltData4 = melt(Data4, id=c('Method'))
    
    
    p1 = ggplot(MeltData1) +
      geom_boxplot(aes(x=Method, y=value, color=Method))+
      facet_grid(.~variable) + 
      theme(axis.text=element_text(size=8),axis.title=element_text(size=10,face="bold")) + 
      ggtitle(paste(Model," graph, p=",Dimension, ", n=",Samplesize[1],sep="")) +
      labs(x=NULL) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(color = "black", size=8, vjust=.8, hjust=0.8),
            axis.ticks.x=element_blank(),
            legend.title=element_text(size=5),
            legend.text=element_text(size=4),
            legend.key.size = unit(1,"line")) +
      scale_x_discrete(labels=NULL)
    
    p2 = ggplot(MeltData2) +
      geom_boxplot(aes(x=Method, y=value, color=Method))+
      facet_grid(.~variable) + 
      theme(axis.text=element_text(size=8),axis.title=element_text(size=10,face="bold")) + 
      ggtitle(paste(Model," graph, p=",Dimension, ", n=",Samplesize[2],sep="")) +
      labs(x=NULL) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(color = "black", size=8, vjust=.8, hjust=0.8),
            axis.ticks.x=element_blank(),
            legend.title=element_text(size=5),
            legend.text=element_text(size=4),
            legend.key.size = unit(1,"line")) +
      scale_x_discrete(labels=NULL)
    
    p3 = ggplot(MeltData3) +
      geom_boxplot(aes(x=Method, y=value, color=Method))+
      facet_grid(.~variable) + 
      theme(axis.text=element_text(size=8),axis.title=element_text(size=10,face="bold")) + 
      ggtitle(paste(Model," graph, p=",Dimension, ", n=",Samplesize[3],sep="")) +
      labs(x=NULL) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(color = "black", size=8, vjust=.8, hjust=0.8),
            axis.ticks.x=element_blank(),
            legend.title=element_text(size=5),
            legend.text=element_text(size=4),
            legend.key.size = unit(1,"line")) +
      scale_x_discrete(labels=NULL)
    
    p4 = ggplot(MeltData4) +
      geom_boxplot(aes(x=Method, y=value, color=Method))+
      facet_grid(.~variable) + 
      theme(axis.text=element_text(size=8),axis.title=element_text(size=10,face="bold")) + 
      ggtitle(paste(Model," graph, p=",Dimension, ", n=",Samplesize[4],sep="")) +
      labs(x=NULL) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(color = "black", size=8, vjust=.8, hjust=0.8),
            axis.ticks.x=element_blank(),
            legend.title=element_text(size=5),
            legend.text=element_text(size=4),
            legend.key.size = unit(1,"line")) +
      scale_x_discrete(labels=NULL)
    
    scale=0.8
    
    pdf(paste(Dimension,Model,Neigh[i],"NeighborhoodDiagnostics.pdf",sep=""), width = 10*scale, height = 7*scale)
    grid.arrange(p1, p2, p3, p4, ncol=2)
    dev.off() 
    
  }
  
}
