
# non-central chi-squared density to be used for optimization
chisq.dens <- function(ncp, z, df){
  res <- dchisq(x=z, df=df, ncp=ncp, log=TRUE)
  return(res)
}

# density of the F-distribution to be used for optimization
F.dens  <- function(ncp, Fstat, df1, df2){
  res <- df(x=Fstat, df1=df1, df2=df2, ncp=ncp, log=TRUE)
  return(res)
}

# log density of the noncentral hypergeometric distribution to be used for optimization
noncentralhypergeom.log.dens <- function(or, x){
  n11 <- x[1,1]
  n12 <- x[1,2]
  n21 <- x[2,1]
  n22 <- x[2,2]
  # row sums of table
  np1 <- n11 + n21
  np2 <- n12 + n22
  # column sums of table
  n1p <- n11 + n12
  n2p <- n21 + n22
  # sum of all entries in the table
  n <- np1 + np2
  
  # use efficient implementation of non-centr. hypergeom. distr.
  res <- log(dnoncenhypergeom(x=n11, n1=n1p, n2=n2p, m1=np1, psi=or))
  return (res)
}

# log density under two-sided alternative for 2x2 table
log.dens.alternative <- function(or, x){
  log(exp(noncentralhypergeom.log.dens(or=or, x=x)) + exp(noncentralhypergeom.log.dens(or=1/or, x=x)))
}

#####################################################################

# BF_10 here (reverse orientation as for minBF)
# normal lik., BF based on distr. of abs(t), which has a folded normal distr.
bf.fn <- function(mu, p, log=FALSE){
  t.star <- qnorm(p/2, lower.tail = FALSE)
  log.num <- log(dnorm(t.star-mu)+dnorm(t.star+mu))
  log.den <- log(2) + dnorm(t.star, log = TRUE)
  if(log==TRUE)
    return(log.num- log.den)
  else
    return(exp(log.num- log.den))
}


# t-distr. with degrees of freedom df (for small samples)
# BF_10 based on the distr. of abs(t), which has a folded t-distr.
bf.ft <- function(mu, p, df, log=FALSE){
  t.star <- qt(p/2, df=df, lower.tail = FALSE)
  log.num <- log(dt(t.star-mu, df=df)+dt(t.star+mu, df=df))
  log.den <- log(2) + dt(t.star, df=df, log = TRUE)
  if(log==TRUE)
    return(log.num - log.den)
  else
    return(exp(log.num - log.den))
}

####################################################################

# transform to logarithm with different bases, input is the log.minBF
transf <- function(x, fun="id"){
    if(fun=="id")
      result <- exp(x)
    if(fun=="log")
      result <- x
    if(fun=="log2"){
      c <- (log(2))^{-1}
      result <- c*x
    }
    if(fun=="log10"){
      c <- (log(10))^{-1}
      result <- c*x
    }
    return(result)
  # }
}

