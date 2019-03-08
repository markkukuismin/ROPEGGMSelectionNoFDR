
ROPEPoissonEED = function(Y,lambda,alpha,breaks = NULL,target=NULL,cv=T){
  
  # Y = n x p data matrix
  
  # lambda = a sequence of penalization parameters for ROPE. The default is five-fold Crossvalidation
  
  # alpha = The level of significance
  
  # breaks = Number of cells in the histogram
  
  if(is.null(breaks)) breaks = "Sturges" # R's default method
  
  n = nrow(Y)
  p = ncol(Y)
  
  if(is.null(target)) target = diag(0,p)
  
  I = diag(1,p)
  
  S = crossprod(Y)/(n-1)
  
  cv.values = NA
  
  if(cv==T){
    ROPECV = crossvalidationROPE(Y,lambda,target=target)
    lambda = ROPECV$cv.rho
    cv.values = ROPECV$crossvalues
  }
  
  hatTheta = ROPE(S,lambda,target=target)

  R = -Scale(hatTheta)

  rho = R[upper.tri(R)]

  EED = -n*log(1 - rho^2)

  #####################################################################################

  # Fit a Poisson regression curve to the empirical distribution of EED

  EEDHist = hist(EED,plot=F,breaks = breaks)
  y = sort(EED)

  x = EEDHist$mids

  s = EEDHist$counts

  malli = glm(s~x,family=poisson(link = "log"))

  Beta <<- coef(malli)

  f = function(y){ exp(Beta[[1]] + Beta[[2]]*y)} # Poisson model
  
  A = integrate(f,lower=0,upper=Inf)$value
  
  yfit = f(y)/A

  beta0 = Beta[[1]]; beta1 = Beta[[2]]
  
  #####

  B = (1/A)*exp(beta0)/beta1

  # Determine the 100(1-alpha) percentile 

  EED = -n*log(1 - (R*(1-I))^2)

  a = log((1 - alpha)/B + 1)*beta1^(-1)
  
  IndROPE = ifelse(EED < a,0,1)
  
  diag(IndROPE) = 0
  
  ###########################################
  
  
  results = list("IndROPEPoisson" = IndROPE, "cv.values" = cv.values ,"lambda"=lambda, "EED"=EED, 
                 "yfit"=yfit, "Beta"=Beta)
  
  return(results)
   
}