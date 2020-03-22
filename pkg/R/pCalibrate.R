
pCalibrate <- function(p,  type="exploratory", transform="id"){
  if(min(p) <=0 || max(p) >1)
    stop("All elements of p must lie in (0,1]!")
  
  if(type=="exploratory"){
    log.minBF <- ifelse(log(p)< -1, 1+ log(p) + log(-log(p)), 0)
  }
  if(type=="confirmatory"){
    log.minBF <- ifelse(p<1-exp(-1),  1 + log(1-p) + log(-log(1-p)), 0)
  }
  result <- transf(log.minBF, fun=transform)
  return(result)
}






