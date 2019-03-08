
# 28.02.2019

library(camel)
library(Matrix)
library(spatstat)
library(huge)
library(igraph)
library(GGMridge)
library(rags2ridges)
library(FastGGM)

source("ROPE.txt")
source("crossvalidationROPE.txt")
source("ROPEEEDKernelMethod.r")
source("ROPEEEDEmpiricPoisson.r")

source("Scale.r")
source("GGMRidge.txt")
source("SpSeFallPreMCC.txt")

Seed = read.table("SeedNumber.txt",header=T)
Seed = Seed$x

Model = "Hub" # Hub or Cluster

u = 0.1
v = 0.3

p = 200 # 200 or 500

g = 3 # Nmb of hubs

lambda = seq(0.01,10,length.out=40) # The largest value should be even larger when p increases

M = 50 # The number of simulation rounds

alpha = 0.05 # 100*(1-alpha)% confidence level

B = NULL # Number of breaks in the histogram. Default "Sturges"

Samplesize = c(floor(0.1*p),floor(0.25*p),floor(0.5*p),floor(1*p))

for(n in Samplesize){
  
  nmb = p*(p-1)/2
  
  I = diag(1,p)
  
  #####################################################################################
  
  ResultsROPEPoisson = rep(NA,M)
  ResultsROPEKernel = rep(NA,M)
  ResultsGGMRidge = rep(NA,M)
  ResultsRagsRidge = rep(NA,M)
  ResultsTIGERStARS = rep(NA,M)
  ResultsFastGGM = rep(NA,M)
  
  DiagnosticResultsROPEPoisson = matrix(0,M,5)
  DiagnosticResultsROPEKernel = matrix(0,M,5)
  DiagnosticResultsGGMRidge = matrix(0,M,5)
  DiagnosticResultsRagsRidge = matrix(0,M,5)
  DiagnosticResultsTIGERStARS = matrix(0,M,5)
  DiagnosticResultsFastGGM = matrix(0,M,5)
  
  colnames(DiagnosticResultsROPEPoisson) = c("Sp","Sen","Fall","Pre","MCC")
  colnames(DiagnosticResultsROPEKernel) = colnames(DiagnosticResultsROPEPoisson)
  colnames(DiagnosticResultsGGMRidge) = colnames(DiagnosticResultsROPEPoisson)
  colnames(DiagnosticResultsRagsRidge) = colnames(DiagnosticResultsROPEPoisson)
  colnames(DiagnosticResultsTIGERStARS) = colnames(DiagnosticResultsROPEPoisson)
  colnames(DiagnosticResultsFastGGM) = colnames(DiagnosticResultsROPEPoisson)
  
  set.seed(Seed)
  
  if(Model == "Hub"){
    
    L = huge.generator(n=2*n,d=p,g=g,graph="hub",v=v,u=u)  # Generate data with hub structures
    
    Sigma = L$sigma
    Theta = L$theta
    Theta = as.matrix(Theta)
    
  }
  
  if(Model == "Cluster"){
    
    L = huge.generator(n=2*n,d=p,g=g,graph="cluster",v=v,u=u)  # Generate data with cluster structures
    
    Sigma = L$sigma
    Theta = L$theta
    Theta = as.matrix(Theta)
    
  }
  
  CorR = Scale(Sigma)
  
  summary(CorR[upper.tri(CorR)])
  
  TrueA = L$theta
  TrueA = as.matrix(TrueA)
  TrueG = graph.adjacency(TrueA,mode="undirected",diag=F)
  
  ClustersTrue = walktrap.community(TrueG)
  TrueNmbOfClusters = length(table(membership(ClustersTrue)))
  
  if(Model == "Hub"){
    
    DegreeROPEPoisson = matrix(NA,p,M)
    DegreeROPEKernel = matrix(NA,p,M)
    
    DegreeGGMRidge = matrix(NA,p,M)
    DegreeRagsRidge = matrix(NA,p,M)
    DegreeTIGER = matrix(NA,p,M)
    DegreeFastGGM = matrix(NA,p,M)
    
    Specificity= array(0,c(M,6,g))
    Sensitivity = array(0,c(M,6,g))
    Fallout = array(0,c(M,6,g))
    Precision = array(0,c(M,6,g))
    MCC = array(0,c(M,6,g))
    
    colnames(Specificity) = c("ROPEPoisson","ROPEKer","GGMRidge","RagsRidge","TIGER","FastGGM")
    colnames(Sensitivity) = colnames(Specificity)
    colnames(Fallout) = colnames(Specificity)
    colnames(Precision) = colnames(Specificity)
    colnames(MCC) = colnames(Specificity)
    
  }
  
  set.seed(Seed)
  
  for(i in 1:M){
    
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
    
    # How the other precision matrix estimators hold up?
    
    # Method similar to ROPE, "RagsRidges" which utilizes partial correlation and Strimmer method:
    
    OPT = optPenalty.LOOCV(Yvalid, lambdaMin = .01, lambdaMax = 30, step = 40,verbose = F,graph=F,target = I)
    
    ## Determine support regularized (standardized) precision under optimal penalty
    
    ThetaRagsRidge = ridgeP(cov(Y),lambda=OPT$optLambda,target=I)
    
    ThetaRagsRidge = sparsify(ThetaRagsRidge, threshold = "localFDR",verbose = F)
    
    ARagsRidge = ifelse(ThetaRagsRidge$sparsePrecision != 0, 1, 0); diag(ARagsRidge) = 0
    
    # GGMRidge which uses Fisher Z and Efron method:
    
    GGMRSelect = GGMRidge(Yvalid,cv=T)
    
    GGMRidgeLambda = GGMRSelect$lambda
    
    GGMRidgeR = GGMRidge(Y,cv=F,lambda=GGMRidgeLambda)
    
    GGMRidgeR = GGMRidgeR$GGMRidgePartialCor
    
    AGGMRidge = ifelse(GGMRidgeR != 0, 1, 0); diag(AGGMRidge) = 0
    
    # TIGER with StARS:
    
    out.tiger = camel.tiger(Yvalid,method="slasso",verbose=F) # Let tiger choose the optimal tuning parameter
    lambdaTIGER = camel.tiger.select(out.tiger,criterion="stars",stars.thresh = 0.05,verbose=F)$opt.lambda
    
    ThetaTIGER = camel.tiger(Y,method="slasso",verbose=F,lambda=lambdaTIGER)$icov[[1]]
    ThetaTIGER = as.matrix(ThetaTIGER)
    
    ATIGERStARS = ifelse(ThetaTIGER != 0, 1, 0); diag(ATIGERStARS) = 0
    
    # FastGGM:
    
    # No validations available:
    
    FGGM = FastGGM_Parallel(Y)
    
    AFastGGM = ifelse(FGGM$CI_low_parCor*FGGM$CI_high_parCor < 0, 0, 1)
    diag(AFastGGM) = 0
    
    ###################################################
    
    PoissonG = graph.adjacency(AROPEPoisson,mode="undirected",diag=F)
    KernelG = graph.adjacency(AROPEKernel,mode="undirected",diag=F)
    RagsRidgeG = graph.adjacency(ARagsRidge,mode="undirected",diag=F)
    GGMRidgeG = graph.adjacency(AGGMRidge,mode="undirected",diag=F)
    TIGERG = graph.adjacency(ATIGERStARS,mode="undirected",diag=F)
    FastGGMG = graph.adjacency(AFastGGM,mode="undirected",diag=F)
    
    ClustersPoisson = walktrap.community(PoissonG)
    ClustersKernel = walktrap.community(KernelG)
    ClustersGGMRidge = walktrap.community(GGMRidgeG)
    ClustersRagsRidge = walktrap.community(RagsRidgeG)
    ClustersTIGER = walktrap.community(TIGERG)
    ClustersFastGGM = walktrap.community(FastGGMG)
    
    ResultsROPEPoisson[i] = length(table(membership(ClustersPoisson)))
    ResultsROPEKernel[i] = length(table(membership(ClustersKernel)))
    ResultsGGMRidge[i] = length(table(membership(ClustersGGMRidge)))
    ResultsRagsRidge[i] = length(table(membership(ClustersRagsRidge)))
    ResultsTIGERStARS[i] = length(table(membership(ClustersTIGER)))
    ResultsFastGGM[i] = length(table(membership(ClustersFastGGM)))
    
    if(Model=="Cluster"){
      
      DiagnosticResultsROPEPoisson[i,] = unlist(Diagnostic(AROPEPoisson,TrueA))
      DiagnosticResultsROPEKernel[i,] = unlist(Diagnostic(AROPEKernel,TrueA))
      DiagnosticResultsGGMRidge[i,] = unlist(Diagnostic(AGGMRidge,TrueA))
      DiagnosticResultsRagsRidge[i,] = unlist(Diagnostic(ARagsRidge,TrueA))
      DiagnosticResultsTIGERStARS[i,] = unlist(Diagnostic(ATIGERStARS,TrueA))
      DiagnosticResultsFastGGM[i,] = unlist(Diagnostic(AFastGGM,TrueA))
      
    }
    
    if(Model=="Hub"){
      
      DegreeTrue = colSums(TrueA)
      TrueHubs = order(DegreeTrue, decreasing = T)[1:g]
      
      DegreeROPEPoisson[,i] = colSums(AROPEPoisson)
      DegreeROPEKernel[,i] = colSums(AROPEKernel)
      
      DegreeGGMRidge[,i] = colSums(AGGMRidge)
      DegreeRagsRidge[,i] = colSums(ARagsRidge)
      DegreeTIGER[,i] = colSums(ATIGERStARS)
      DegreeFastGGM[,i] = colSums(AFastGGM)
      
      ####################################################################
      
      # How the neighborhood of a hub node is estimated?
      
      HubEstimated = order(DegreeROPEPoisson[,i],decreasing = T)[1:g] 
      FoundHubs = intersect(HubEstimated,TrueHubs)
      
      if(length(FoundHubs) != 0){
        for(j in FoundHubs){
          DRes = Diagnostic(AROPEPoisson[-j,j],TrueA[-j,j])
          ind = which(FoundHubs == j)
          Specificity[i,1,ind] = DRes$Sp
          Sensitivity[i,1,ind] = DRes$Sen
          Fallout[i,1,ind] = DRes$Fall
          Precision[i,1,ind] = DRes$Pre
          MCC[i,1,ind] = DRes$MCC
        } 
      }
      
      HubEstimated = order(DegreeROPEKernel[,i],decreasing = T)[1:g] 
      FoundHubs = intersect(HubEstimated,TrueHubs)
      
      if(length(FoundHubs) != 0){
        for(j in FoundHubs){
          DRes = Diagnostic(AROPEKernel[-j,j],TrueA[-j,j])
          ind = which(FoundHubs == j)
          Specificity[i,2,ind] = DRes$Sp
          Sensitivity[i,2,ind] = DRes$Sen
          Fallout[i,2,ind] = DRes$Fall
          Precision[i,2,ind] = DRes$Pre
          MCC[i,2,ind] = DRes$MCC
        } 
      }
      
      HubEstimated = order(DegreeGGMRidge[,i],decreasing = T)[1:g] 
      FoundHubs = intersect(HubEstimated,TrueHubs)
      
      if(length(FoundHubs) != 0){
        for(j in FoundHubs){
          DRes = Diagnostic(AGGMRidge[-j,j],TrueA[-j,j])
          ind = which(FoundHubs == j)
          Specificity[i,3,ind] = DRes$Sp
          Sensitivity[i,3,ind] = DRes$Sen
          Fallout[i,3,ind] = DRes$Fall
          Precision[i,3,ind] = DRes$Pre
          MCC[i,3,ind] = DRes$MCC
        } 
      }
      
      HubEstimated = order(DegreeRagsRidge[,i],decreasing = T)[1:g] 
      FoundHubs = intersect(HubEstimated,TrueHubs)
      
      if(length(FoundHubs) != 0){
        for(j in FoundHubs){
          DRes = Diagnostic(ARagsRidge[-j,j],TrueA[-j,j])
          ind = which(FoundHubs == j)
          Specificity[i,4,ind] = DRes$Sp
          Sensitivity[i,4,ind] = DRes$Sen
          Fallout[i,4,ind] = DRes$Fall
          Precision[i,4,ind] = DRes$Pre
          MCC[i,4,ind] = DRes$MCC
        } 
      }
      
      HubEstimated = order(DegreeTIGER[,i],decreasing = T)[1:g] 
      FoundHubs = intersect(HubEstimated,TrueHubs)
      
      if(length(FoundHubs) != 0){
        for(j in FoundHubs){
          DRes = Diagnostic(ATIGERStARS[-j,j],TrueA[-j,j])
          ind = which(FoundHubs == j)
          Specificity[i,5,ind] = DRes$Sp
          Sensitivity[i,5,ind] = DRes$Sen
          Fallout[i,5,ind] = DRes$Fall
          Precision[i,5,ind] = DRes$Pre
          MCC[i,5,ind] = DRes$MCC
        } 
      }
      
      HubEstimated = order(DegreeFastGGM[,i],decreasing = T)[1:g] 
      FoundHubs = intersect(HubEstimated,TrueHubs)
      
      if(length(FoundHubs) != 0){
        for(j in FoundHubs){
          DRes = Diagnostic(AFastGGM[-j,j],TrueA[-j,j])
          ind = which(FoundHubs == j)
          Specificity[i,6,ind] = DRes$Sp
          Sensitivity[i,6,ind] = DRes$Sen
          Fallout[i,6,ind] = DRes$Fall
          Precision[i,6,ind] = DRes$Pre
          MCC[i,6,ind] = DRes$MCC
        } 
      }
      
      
    }
    
    cat("\r", 100*round(i/M,2), "% complete")
    
  }
  
  
  ##############################################################################
  
  # What about hub-nodes, how many were identified right?
  
  # Now there are g hubs:
  
  # True hubs: 
  
  if(Model == "Hub"){
    
    DegreeTrue = colSums(TrueA)
    TrueHubs = order(DegreeTrue, decreasing = T)[1:g]
    
    NmbOfFoundHubs = matrix(NA,M,6)
    colnames(NmbOfFoundHubs) = c("ROPEPoisson","ROPEKernel","GGMRidge","RagsRidge","TIGER","FastGGM")
    
    for(i in 1:M){
      
      if(sum(DegreeROPEPoisson[,i])!=0){
        HubEstimated = order(DegreeROPEPoisson[,i],decreasing = T)[1:g] 
        FoundHubs = intersect(HubEstimated,TrueHubs)
        NmbOfFoundHubs[i,1] = length(FoundHubs)
      }else{
        NmbOfFoundHubs[i,1] = 0
      }
      
      if(sum(DegreeROPEKernel[,i])!=0){
        HubEstimated = order(DegreeROPEKernel[,i],decreasing = T)[1:g] 
        FoundHubs = intersect(HubEstimated,TrueHubs)
        NmbOfFoundHubs[i,2] = length(FoundHubs)
      }else{
        NmbOfFoundHubs[i,2] = 0
      }
      
      
      if(sum(DegreeGGMRidge[,i])!=0){
        HubEstimated = order(DegreeGGMRidge[,i],decreasing = T)[1:g] 
        FoundHubs = intersect(HubEstimated,TrueHubs)
        NmbOfFoundHubs[i,3] = length(FoundHubs) 
      } else{
        NmbOfFoundHubs[i,3] = 0
      }
      
      if(sum(DegreeRagsRidge[,i])!=0){
        HubEstimated = order(DegreeRagsRidge[,i],decreasing = T)[1:g] 
        FoundHubs = intersect(HubEstimated,TrueHubs)
        NmbOfFoundHubs[i,4] = length(FoundHubs) 
      }else{
        NmbOfFoundHubs[i,4] = 0
      }
      
      if(sum(DegreeTIGER[,i])!=0){
        HubEstimated = order(DegreeTIGER[,i],decreasing = T)[1:g] 
        FoundHubs = intersect(HubEstimated,TrueHubs)
        NmbOfFoundHubs[i,5] = length(FoundHubs) 
      }else{
        NmbOfFoundHubs[i,5] = 0
      }
      
      if(sum(DegreeFastGGM[,i])!=0){
        HubEstimated = order(DegreeFastGGM[,i],decreasing = T)[1:g] 
        FoundHubs = intersect(HubEstimated,TrueHubs)
        NmbOfFoundHubs[i,6] = length(FoundHubs) 
      }else{
        NmbOfFoundHubs[i,6] = 0
      }
      
    }
    
    GGMMethod = c(rep("ROPEPoisson",times=M),rep("ROPEKernel",times=M),rep("GGMRidge",times=M),
                  rep("RagsRidge",times=M),rep("TIGER",times=M),rep("FastGGM",times=M))
    
    
  }
  
  data = data.frame("ROPEPoisson" = ResultsROPEPoisson,"ROPEKernel"=ResultsROPEKernel,
                    "GGMRidge"=ResultsGGMRidge,"RagsRidge"=ResultsRagsRidge,"TIGERStARS" = ResultsTIGERStARS, 
                    "FastGGM" = ResultsFastGGM,"TrueNmbOfClusters" = rep(TrueNmbOfClusters,M))
  
  #################################################################################################
  
  # These are for the whole graph - in particular, for the cluster graph:
  
  DiagnosticResults = rbind(DiagnosticResultsROPEPoisson,DiagnosticResultsROPEKernel,DiagnosticResultsGGMRidge,
                            DiagnosticResultsRagsRidge,DiagnosticResultsTIGERStARS,DiagnosticResultsFastGGM)
  
  DiagnosticResults = as.data.frame(DiagnosticResults)
  
  DiagnosticResults = cbind(DiagnosticResults,"Method"=rep(c("ROPEPoisson","ROPEKernel","GGMRidge","RagsRidge",
                                                             "TIGER","FastGGM"),each=M))
  
  #################################################################################################
  
  if(Model=="Hub"){
    
    FirstHubDiagnostics = data.frame(Sp=rep(0,6*M),Sen=rep(0,6*M),Fall=rep(0,6*M),Pre=rep(0,6*M),
                                     MCC=rep(0,6*M),
                                     Method=rep(c("ROPEPoisson","ROPEKernel","GGMRidge","RagsRidge","TIGER",
                                                  "FastGGM"),each=M)
    )
    
    SecondHubDiagnostics = FirstHubDiagnostics
    ThirdHubDiagnostics = FirstHubDiagnostics
    
    FirstHubDiagnostics$Sp  = as.vector(Specificity[,1:6,1])
    FirstHubDiagnostics$Sen = as.vector(Sensitivity[,1:6,1])
    FirstHubDiagnostics$Fall= as.vector(Fallout[,1:6,1])
    FirstHubDiagnostics$Pre = as.vector(Precision[,1:6,1])
    FirstHubDiagnostics$MCC = as.vector(MCC[,1:6,1])
    
    SecondHubDiagnostics$Sp  = as.vector(Specificity[,1:6,2])
    SecondHubDiagnostics$Sen = as.vector(Sensitivity[,1:6,2])
    SecondHubDiagnostics$Fall= as.vector(Fallout[,1:6,2])
    SecondHubDiagnostics$Pre = as.vector(Precision[,1:6,2])
    SecondHubDiagnostics$MCC = as.vector(MCC[,1:6,2])
    
    ThirdHubDiagnostics$Sp  = as.vector(Specificity[,1:6,3])
    ThirdHubDiagnostics$Sen = as.vector(Sensitivity[,1:6,3])
    ThirdHubDiagnostics$Fall= as.vector(Fallout[,1:6,3])
    ThirdHubDiagnostics$Pre = as.vector(Precision[,1:6,3])
    ThirdHubDiagnostics$MCC = as.vector(MCC[,1:6,3])
    
  }
  
  
  if(Model == "Cluster"){
    
    write.csv(DiagnosticResults,file=paste(Model,"WholeGraphDiagnosticForDimension",p,"n=",n, ".csv", sep=""),
              row.names=F) 
    
  }
  
  write.csv(data,file=paste(Model,"NmbResultsForDimension",p,"n=",n, ".csv", sep=""),row.names=F)
  
  if(Model == "Hub"){
    
    write.csv(NmbOfFoundHubs,file=paste(Model,"NmbOfCorretlyIdentHubsForDimension",p,"n=",n,".csv",sep=""),row.names=F)
    
    write.csv(FirstHubDiagnostics,file=paste(Model,"FirstNeighborhoodDiagnosticsForDimension",p,"n=",n,".csv",sep=""),
              row.names=T)
    
    write.csv(SecondHubDiagnostics,file=paste(Model,"SecondNeighborhoodDiagnosticsForDimension",p,"n=",n,".csv",sep=""),
              row.names=T)
    
    write.csv(ThirdHubDiagnostics,file=paste(Model,"ThirdNeighborhoodDiagnosticsForDimension",p,"n=",n,".csv",sep=""),
              row.names=T)
    
  }
  
  
}

par(mfrow=c(2,3))

f = function(m) t(m)[nrow(m):1,] # For plotting the images
plot(im(f(Theta)), main="Real Adjacency",col = gray(seq(1,0,length.out=10)))
plot(im(f(AROPEPoisson)), main="ROPE with EED (Empiric Poisson)",col = gray(seq(1,0,length.out=10)))
plot(im(f(AROPEKernel)), main="ROPE with EED (Kernel)",col = gray(seq(1,0,length.out=10)))
plot(im(f(ARagsRidge)), main="RagsRidge",col = gray(seq(1,0,length.out=10)))
plot(im(f(ATIGERStARS)), main="TIGER adjacency",col = gray(seq(1,0,length.out=10)))
plot(im(f(AFastGGM)), main="FastGGM adjacency",col = gray(seq(1,0,length.out=10)))

# How many percentages of estimated clusters equals to the "true" number of clusters?

round(100*sum(ResultsROPEPoisson == TrueNmbOfClusters)/M,2)
round(100*sum(ResultsROPEKernel == TrueNmbOfClusters)/M,2)
round(100*sum(ResultsGGMRidge == TrueNmbOfClusters)/M,2)
round(100*sum(ResultsRagsRidge == TrueNmbOfClusters)/M,2)
round(100*sum(ResultsTIGERStARS == TrueNmbOfClusters)/M,2)
round(100*sum(ResultsFastGGM == TrueNmbOfClusters)/M,2)

par(mfrow=c(1,1))
boxplot(ResultsROPEPoisson,ResultsROPEKernel,ResultsGGMRidge,ResultsRagsRidge,ResultsTIGERStARS,ResultsFastGGM,
        names = c("Poisson","Kernel","GGMRidge","RagsRidge","TIGERStARS","FastGGM"),
        ylim=c(0,max(c(ResultsTIGERStARS,ResultsGGMRidge,ResultsROPEKernel,ResultsROPEPoisson,ResultsRagsRidge,
                       ResultsFastGGM))),
        main="Nmb of estimated clusters")
abline(h=TrueNmbOfClusters,col="red",lty=2,lwd=2)

ExtraNmb = as.vector(NmbOfFoundHubs)

DTable = table(ExtraNmb,GGMMethod)

barplot(DTable, main="Hub-node distribution by Method and nmb of correctly identified hubs",
        xlab="Method", col=rainbow(length(ExtraNmb)),ylab="Counts",
        legend = rownames(DTable), beside=TRUE) 

summary(DiagnosticResultsROPEPoisson)
summary(DiagnosticResultsROPEKernel)

summary(DiagnosticResultsGGMRidge)
summary(DiagnosticResultsRagsRidge)
summary(DiagnosticResultsTIGERStARS)
summary(DiagnosticResultsFastGGM)