Scale= function(A){
  
  p = ncol(A)
  
  diagA = diag(diag(A)^(-.5),p)
  
  outer(diag(diagA), diag(diagA)) * A
  
  
}