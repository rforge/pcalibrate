

# transformation of BF to posterior probability (of the hypothesis in the numerator of the bf)
BF2pp <- function(BF, prior.prob=0.5){
  prior.odds <- prior.prob/(1-prior.prob)
  log.post.odds <- log(prior.odds) + log(BF)
  post.prob <- exp(log.post.odds)/(1+ exp(log.post.odds))
  return (post.prob)
}

           
########################################
# p-based Bayes factors
########################################


pCalibrate <- function(p,  alternative="noninformative", transform="id"){
  if(min(p) <=0 || max(p) >1)
    stop("All elements of p must lie in (0,1]!")
  
  if(alternative=="noninformative"){
    log.minBF <- ifelse(log(p)< -1, 1+ log(p) + log(-log(p)), 0)
  }
  if(alternative=="informative"){
    log.minBF <- ifelse(p<1-exp(-1),  1 + log(1-p) + log(-log(1-p)), 0)
  }
  result <- transf(log.minBF, fun=transform)
  return(result)
}

   
         
###############################################
# Bayes factors based on a test statistic
###############################################


# calibrations using the z-value
zCalibrate  <-  function(p, type="two.sided", alternative="normal", transform="id"){
  if(min(p) <=0 || max(p) >1)
    stop("All elements of p must lie in (0,1]!")
  
  if(type=="one.sided" && alternative=="simple"){
    t <- qnorm(p, lower.tail=FALSE)
    log.minBF <- ifelse(t<=0, 0, -0.5*t^2)
  }
  if(type=="two.sided" && alternative=="simple"){
    t.star <- qnorm(p/2, lower.tail=FALSE)
    log.minBF <- numeric(length(p))
    for(i in 1:length(p)){
      res <- optimize(bf.fn, p=p[i], log=TRUE, lower=0, upper=t.star[i]+1, maximum=TRUE)
      log.minBF[i] <- - res$objective
    }
  }
  if(alternative=="normal"){
    t <- qnorm(p/2, lower.tail=FALSE)
      log.minBF <- ifelse(t<=1, 0, 1/2 + log(t) - t^2/2)
  }
  if(alternative=="local"){
    target <- function(K, t){
      left <- K*(dnorm(K+t)+dnorm(K-t))
      right <- pnorm(K-t) - pnorm(-(K+t))
      return((left-right)^2)
    }
    
    t <- qnorm(p/2, lower.tail=FALSE)
    log.minBF <- rep(NA, length(p))
    
    for(i in 1:length(p)){
      if(t[i] <= 1) log.minBF[i] <- 0
      if(t[i] > 1){
         K <- optimize(target, lower=0, upper=4*t[i], t=t[i])$minimum
         log.minBF[i] <- log(2) + dnorm(t[i], log=TRUE) - log(dnorm(K+t[i])+dnorm(K-t[i]))
      }
    }
  }
  result <- transf(log.minBF, fun=transform)
  return(result)
}


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
  }
  result <- transf(log.minBF, fun=transform)
  return(result)
}



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
  }
  if(alternative=="normal"){
      log.minBF <- FCalibrate(p=p, n=n, d=1, alternative="chi.squared", transform="log")
  }
  result <- transf(log.minBF, fun=transform)
  return(result)
}


LRCalibrate <- function(p, df, alternative="gamma", transform="id"){
  if(min(p) <=0 || max(p) >1)
    stop("All elements of p must lie in (0,1]!")
  
  if(alternative=="simple"){
    k <- length(p)
    m <- length(df)
    log.minBF <- matrix(NA, ncol=k, nrow=m)
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
    k <- length(p)
    m <- length(df)
    log.minBF <- matrix(NA, ncol=k, nrow=m)
    for(i in 1:k){
      for(j in 1:m){
        z <- qchisq(p[i], df=df[j], lower.tail=FALSE)
        log.minBF[j, i] <- ifelse(z>= df[j], df[j]/2*log(z/df[j])-(z-df[j])/2, 0)
      }
    }
  }
  result <- transf(log.minBF, fun=transform)
  return(result)
}

#########################################################

#############################################
# minBFs for 2x2 contingency tables
#############################################


# 2x2 contingency table should be a matrix
twoby2Calibrate <- function(x, type="two.sided", alternative="normal", direction=NULL, transform.bf="id"){
  
  min.entry <- min(x)
  if(min.entry < 0){
    stop("All entries of x must be non-negative.")
  }
  
  if(type=="one.sided" && direction=="greater"){
    p.fi <- fisher.test(x=x, or=1, alternative = "greater")$p.value
    
    a <- x[1,1]
    b <- x[1,2]
    c <- x[2,1]
    d <- x[2,2]
    p.mid <- p.fi - dhyper(x=a, m=a+b, n=c+d, k=a+c)/2
    
    I2 <- matrix(data=c(1, 0, 0, 1), nrow=2, byrow=TRUE)
    tab <- x + I2
    p.lie <- exact2x2(x=tab, or = 1, alternative = "greater",
                        conf.int = TRUE, conf.level = 0.95,
                        conditional=TRUE, paired=FALSE)$p.value
    p.value <- c(p.fi=p.fi, p.mid=p.mid, p.lie=p.lie)
    
    if(min.entry ==0){
      log.minBF <- NA
      warning("The minimum Bayes factor is not available for a table with entries =0!")
    }
    if(min.entry >0){
    log.num <- noncentralhypergeom.log.dens(x=x, or=1) 
    n <- x[1,1] + x[1,2] + x[2,1] + x[2,2]
    upper <- (n/2 -1)^2
    log.denom <- optimize(noncentralhypergeom.log.dens, x=x, lower=1, 
                          upper=upper, maximum=TRUE)$objective
    log.minBF <- min(log.num-log.denom, 0)
    }
  }
  if(type=="one.sided" && direction=="less"){
    p.fi <- fisher.test(x=x, or=1, alternative = "less")$p.value
    
    a <- x[1,1]
    b <- x[1,2]
    c <- x[2,1]
    d <- x[2,2]
    p.mid <- p.fi - dhyper(x=a, m=a+b, n=c+d, k=a+c)/2
    
    I2 <- matrix(data=c(1, 0, 0, 1), nrow=2, byrow=TRUE)
    tab <- x + I2
    p.lie <- 1 - exact2x2(x=tab, or = 1, alternative = "greater",
                      conf.int = TRUE, conf.level = 0.95,
                      conditional=TRUE, paired=FALSE)$p.value
    p.value <- c(p.fi=p.fi, p.mid=p.mid, p.lie=p.lie)
    
    if(min.entry ==0){
      log.minBF <- NA
      warning("The minimum Bayes factor is not available for a table with entries =0!")
    }
    if(min.entry >0){
    log.num <- noncentralhypergeom.log.dens(x=x, or=1) 
    n <- x[1,1] + x[1,2] + x[2,1] + x[2,2]
    upper <- (n/2 -1)^2
    log.denom <- optimize(noncentralhypergeom.log.dens, x=x, lower=1/upper, 
                          upper=1, maximum=TRUE)$objective
    theta.hat <- optimize(noncentralhypergeom.log.dens, x=x, lower=1/upper, 
                          upper=1, maximum=TRUE)$maximum
    log.minBF <- min(log.num-log.denom, 0)
    }
  }
  if(type=="two.sided"){
  # probability-based method
  p.pb <- exact2x2(x=x, or = 1, alternative = "two.sided",
                    tsmethod = "minlike", conf.int = TRUE, conf.level = 0.95,
                    conditional=TRUE, paired=FALSE)$p.value
  # central P-value
  p.ce <- exact2x2(x=x, or = 1, alternative = "two.sided",
                    tsmethod = "central", conf.int = TRUE, conf.level = 0.95,
                    conditional=TRUE, paired=FALSE)$p.value
  # Blaker's P-value
  p.bl <- exact2x2(x=x, or = 1, alternative = "two.sided",
                    tsmethod = "blaker", conf.int = TRUE, conf.level = 0.95,
                    conditional=TRUE, paired=FALSE)$p.value
  # mid P-value
  p.mid <- suppressWarnings(tab2by2.test(x, y = NULL,
                        correction = F)$p.value[2,"midp.exact"])
  # Liebermeister's significance measure
  I2 <- matrix(data=c(1, 0, 0, 1), nrow=2, byrow=TRUE)
  tab <- x + I2
  p.minus <- exact2x2(x=tab, or = 1, alternative = "greater",
                      conf.int = TRUE, conf.level = 0.95,
                      conditional=TRUE, paired=FALSE)$p.value
  p.lie <- 2*min(p.minus, 1-p.minus)
  
  p.value <- c(p.pb=p.pb, p.ce=p.ce, p.bl=p.bl, p.mid=p.mid, p.lie=p.lie)
  }
  if(type=="two.sided" && alternative=="simple"){
    if(min.entry ==0){
      log.minBF <- NA
      warning("The minimum Bayes factor is not available for a table with entries =0!")
    }
    if(min.entry >0){
    log.num <- noncentralhypergeom.log.dens(x=x, or=1) + log(2)
    n <- x[1,1] + x[1,2] + x[2,1] + x[2,2]
    upper <- (n/2 -1)^2
    log.denom <- optimize(log.dens.alternative, x=x, lower=1, upper=upper, maximum=TRUE)$objective
    log.minBF <- min(log.num-log.denom, 0)
    }
  }
  if(type=="two.sided" && alternative=="normal"){
    # Li & Clyde (2016) approach using a generalized g-prior
    a <- x[1,1]
    b <- x[1,2]
    c <- x[2,1]
    d <- x[2,2]
    
    if(min.entry ==0){
      log.minBF <- NA
      warning("Minimum Bayes factor is not defined for a table with entries =0!")
    }
    if(min.entry >0){
    # row sums of table
    n1 <- a+b
    n0 <- c+d
    # column sums of table
    m1 <- a+c
    m0 <- b+d
    # sum of all entries in the table
    n <- a+b+c+d
    
    C <- a*b/n1/(a*b/n1 + c*d/n0)
    
    # squared Wald statistic
    Qm <- (log(a)-log(b)-log(c)+log(d))^2 *((1-C)^2*a*b/n1 + C^2*c*d/n0)
    
    # deviance
    zm <- 2*(a*log(a)+b*log(b)+c*log(c)+d*log(d) -m1*log(m1) - m0*log(m0) - n1*log(n1) - n0*log(n0) + n*log(n))
    
    # empirical Bayes estimate for g
    g.hat <- max(Qm-1,0)
    
    if(g.hat==0)
    {
      log.minBF <- min(-zm/2+ 1/2*(log(a*b/n1 + c*d/n0) - log(m1*m0/n)) +
                         1/2*(((1-C)^2*a*b/n1 + C^2*c*d/n0)*(log(a)-log(b)-log(c)+log(d))^2),0)
    } else {
      minBF <- min(((a*b/n1 +c*d/n0)/(m1*m0/n))^{1/2}*((1-C)^2*a*b/n1 + C^2*c*d/n0)^{1/2}*
                     abs(log(a)-log(b)-log(c)+log(d))*exp((1-zm)/2), 1)
      log.minBF <- log(minBF)
      }
     }
  }
  res.bf <- transf(log.minBF, fun=transform.bf)
  return(list(minBF=res.bf, p.value=p.value))
}


#################################################################   

# formatting of Bayes factors

roundBF <- function(BF, digits="default"){
  if(digits=="default"){
    if(BF < 1/1000)
      result <- "< 1/1000"
    if((BF >= 1/1000) & (BF <= 1/10))
      result <- paste("1/", as.character(round(1/BF)), sep="")
    if((BF > 1/10) & (BF < 1))
      result <- paste("1/", as.character(round(1/BF, digits=1)), sep="")
    if((BF < 10) & (BF >= 1))
      result <- as.character(round(BF, digits=1))
    if((BF >= 10) & (BF <= 1000))
      result <- as.character(round(BF))
    if(BF > 1000)
      result <- "> 1000"
  }
  else{
    if(BF < 1)
      result <- paste("1/", as.character(round(1/BF, digits=digits)), sep="")
    else
      result <- as.character(round(BF, digits=digits))
  }
  return(result)
}


formatBF <- function(BF, digits="default"){
  result <- character()
  for(i in 1:length(BF))
    result[i] <- roundBF(BF[i], digits=digits)
  return(result)
}






