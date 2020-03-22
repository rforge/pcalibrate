
# for the linear model with d covariates
# p is a vector, 
# either n or d can be a vector, but the other one must be a scalar
FCalibrate <- function(p, n, d, alternative="chi.squared", intercept=TRUE, transform="id"){
  if(min(p) <=0 || max(p) >1)
    stop("All elements of p must lie in (0,1]!")
  
  if(length(n)>1 && length(d)>1){
    stop("Either n or d must be a scalar.")
  }
  if(intercept==TRUE){
    n1 <- n-1
  }
  if(intercept==FALSE){
    n1 <- n
  }
  if(length(n) > 1){
    k <- length(p)
    m <- length(n)
    if(alternative=="simple"){
      log.minBF <- matrix(NA, ncol=k, nrow=m)
      for(j in 1:m){
        Fstat <- qf(p=p, df1=d, df2=n1[j]-d, lower.tail = FALSE)
        for(i in 1:k){
          if(p[i]==1)
            log.minBF[j,i] <- 0
          else{
            K <- optimize(F.dens, Fstat=Fstat[i], df1=d, df2=n1[j]-d, lower=0, upper=d*Fstat[i]+1,
                          maximum=TRUE)$maximum
            log.minBF[j,i] <- min( F.dens(ncp=0, Fstat=Fstat[i], df1=d, df2=n1[j]-d)- 
                                     F.dens(ncp=K, Fstat=Fstat[i], df1=d, df2=n1[j]-d), 0)
          }
        }
      }
      # }
    }
    if(alternative=="chi.squared"){
      log.minBF <- matrix(NA, ncol=k, nrow=m)
      for(j in 1:m){
        Fstat <- qf(p, df1=d, df2=n1[j]-d, lower.tail = FALSE)
        R2 <- 1/(1+(n1[j]-d)/(d*Fstat))
        log.minBF[j, ] <- ifelse(R2>d/n1[j], n1[j]/2*log(n1[j])+d/2*log(R2/d) 
                                 + (n1[j]-d)/2*(log((1-R2)/(n1[j]-d))), 0)
      }
    }
    if(k==1)
      log.minBF <- as.vector(log.minBF)
  }
  else{
    k <- length(p)
    m <- length(d)
    if(alternative=="simple"){
      log.minBF <- matrix(NA, ncol=k, nrow=m)
      for(j in 1:m){
        Fstat <- qf(p=p, df1=d[j], df2=n1-d[j], lower.tail=FALSE)
        for(i in 1:k){
          if(p[i]==1)
            log.minBF[j,i] <- 0
          else{
            K <- optimize(F.dens, Fstat=Fstat[i], df1=d[j], df2=n1-d[j], 
                          lower=0, upper=d[j]*Fstat[i]+1, maximum=TRUE)$maximum
            log.minBF[j,i] <- min( F.dens(ncp=0, Fstat=Fstat[i], df1=d[j], df2=n1-d[j])- 
                                     F.dens(ncp=K, Fstat=Fstat[i], df1=d[j], df2=n1-d[j]), 0)
          }
        }
      }
    }
    if(alternative=="chi.squared"){
      log.minBF <- matrix(NA, ncol=k, nrow=m)
      for(j in 1:m){
        Fstat <- qf(p, df1=d[j], df2=n1-d[j], lower.tail=FALSE)
        R2 <- 1/(1+(n1-d[j])/(d[j]*Fstat))
        log.minBF[j, ] <- ifelse(R2>d[j]/n1, n1/2*log(n1)+d[j]/2*log(R2/d[j]) 
                                 + (n1-d[j])/2*(log((1-R2)/(n1-d[j]))), 0)
      }
    }
    if(k==1 | m==1)
      log.minBF <- as.vector(log.minBF)
  }
  result <- transf(log.minBF, fun=transform)
  return(result)
}