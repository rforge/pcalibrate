
tCalibrate <- function(p, n, type="two.sided", alternative="normal", transform="id"){
  if(min(p) <=0 || max(p) >1)
    stop("All elements of p must lie in (0,1]!")
  
  if(type=="one.sided" && alternative=="simple"){
    k <- length(p)
    m <- length(n)
    log.minBF <- matrix(NA, ncol=k, nrow=m)
    for(j in 1:m){
      t <- qt(p, df=n[j]-2, lower.tail = FALSE)
      log.minBF[j, ] <- ifelse(t <= 0, 0, -(n[j]-1)/2*log(1+1/(n[j]-2)*t^2))
    }
    if(k==1 | m==1)
      log.minBF <- as.vector(log.minBF)
  }
  if(type=="two.sided" && alternative=="simple"){
    k <- length(p)
    m <- length(n)
    log.minBF <- matrix(NA, ncol=k, nrow=m)
    df <- n-2
    for(j in 1:m){
      for(i in 1:k){
        t.star <- qt(p[i]/2, df=n[j]-2, lower.tail = FALSE)
        res <- optimize(bf.ft, p=p[i], df=df[j], log=TRUE, lower=0, upper=t.star+1, maximum=TRUE)
        log.minBF[j,i] <- - res$objective
      }
    }
    if(k==1 | m==1)
      log.minBF <- as.vector(log.minBF)
  }
  if(alternative=="normal"){
    log.minBF <- FCalibrate(p=p, n=n, d=1, alternative="chi.squared", transform="log")
  }
  
  result <- transf(log.minBF, fun=transform)
  return(result)
}