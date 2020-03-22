
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
    # twice the min. one-sided mid p-value
    # case alternative = "greater"
    p.fi.g <- fisher.test(x=x, or=1, alternative = "greater")$p.value
    a <- x[1,1]
    b <- x[1,2]
    c <- x[2,1]
    d <- x[2,2]
    p.mid.g <- p.fi.g - dhyper(x=a, m=a+b, n=c+d, k=a+c)/2
    
    # case alternative = "less"
    p.fi.l <- fisher.test(x=x, or=1, alternative = "less")$p.value
    p.mid.l <- p.fi.l - dhyper(x=a, m=a+b, n=c+d, k=a+c)/2
    
    p.mid <- 2*min(p.mid.g, p.mid.l)
    
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