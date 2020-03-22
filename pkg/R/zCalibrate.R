
# calibrations using the z-value
zCalibrate  <-  function(p=NULL, z=p2z(p), type="two.sided", alternative="normal", 
                         transform="id"){
  
  if(type=="one.sided" && alternative=="simple"){
    # t <- qnorm(p, lower.tail=FALSE)
    log.minBF <- ifelse(z<=0, 0, -0.5*z^2)
  }
  if(type=="two.sided" && alternative=="simple"){
    # t.star <- qnorm(p/2, lower.tail=FALSE)
    if(is.null(p))
      p <- z2p(z=z, type=type)
    log.minBF <- numeric(length(z))
    for(i in 1:length(z)){
      res <- optimize(bf.fn, p=p[i], log=TRUE, lower=0, upper=z[i]+1, maximum=TRUE)
      log.minBF[i] <- - res$objective
    }
  }
  if(alternative=="normal"){
    # t <- qnorm(p/2, lower.tail=FALSE)
    log.minBF <- ifelse(z<=1, 0, 1/2 + log(z) - z^2/2)
  }
  if(alternative=="local"){
    target <- function(K, t){
      left <- K*(dnorm(K+t)+dnorm(K-t))
      right <- pnorm(K-t) - pnorm(-(K+t))
      return((left-right)^2)
    }
    
    # t <- qnorm(p/2, lower.tail=FALSE)
    log.minBF <- rep(NA, length(z))
    
    for(i in 1:length(z)){
      if(z[i] <= 1) log.minBF[i] <- 0
      if(z[i] > 1){
        K <- optimize(target, lower=0, upper=4*z[i], t=z[i])$minimum
        log.minBF[i] <- log(2) + dnorm(z[i], log=TRUE) - log(dnorm(K+z[i])+dnorm(K-z[i]))
      }
    }
  }
  result <- transf(log.minBF, fun=transform)
  return(result)
}
