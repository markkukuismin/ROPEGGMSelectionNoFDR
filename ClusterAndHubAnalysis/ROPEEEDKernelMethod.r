
ROPEKernelEED = function(Y,lambda,alpha,target=NULL,cv=T){

	library(spatstat)
  
  # Y = n x p data matrix
  
  # lambda = a sequence of penalization parameters for ROPE. The default is five-fold Crossvalidation
  
  # alpha = The level of significance
  
  
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

  # Determine the 100(1-alpha) percentile with kernel density from empiric EED values:
  
  EmpDist = density(EED)
  
  Threshold = quantile(EmpDist,probs=1-alpha)

  EEDMat = -n*log(1 - (R*(1-I))^2)

  # No FDR control
  
  IndROPE = ifelse(EEDMat < Threshold, 0, 1)
  
  diag(IndROPE) = 0
  
  ###########################################
  
  # Calculate P-values, just for fun
  
  EmpiricCDF = CDF(EmpDist)

  Pvalues = EmpiricCDF(EED)
  
  ################
  
  results = list("IndROPEKernel" = IndROPE, "cv.values" = cv.values ,"lambda"=lambda, 
                 "EED"=EED,Pvalues=Pvalues)
  
  return(results)
   
}