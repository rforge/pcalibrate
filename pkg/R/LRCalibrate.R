
LRCalibrate <- function(p, df, alternative="gamma", transform="id"){
  if(min(p) <=0 || max(p) >1)
    stop("All elements of p must lie in (0,1]!")
  
  k <- length(p)
  m <- length(df)
  log.minBF <- matrix(NA, ncol=k, nrow=m)
  
  if(alternative=="simple"){
    for(j in 1:m){
      z <- qchisq(p=p, df=df[j], lower.tail=FALSE)
      for(i in 1:k){
        if(p[i]==1)
          log.minBF[j,i] <- 0
        else{
          K <- optimize(chisq.dens, z=z[i], df=df[j], lower=0, upper=z[i]+1,  maximum=TRUE)$maximum
          log.minBF[j,i] <- min( chisq.dens(ncp=0, z[i], df=df[j])-chisq.dens(ncp=K, z[i], df=df[j]), 0)
        }
      }
    }
  }
  if(alternative=="gamma"){
    for(i in 1:k){
      for(j in 1:m){
        z <- qchisq(p[i], df=df[j], lower.tail=FALSE)
        log.minBF[j, i] <- ifelse(z>= df[j], df[j]/2*log(z/df[j])-(z-df[j])/2, 0)
      }
    }
  }
  
  if(k==1 | m==1)
    log.minBF <- as.vector(log.minBF)
  result <- transf(log.minBF, fun=transform)
  return(result)
}