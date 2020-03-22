
p2z <- function(p, type="two.sided"){
  if(min(p) <=0 || max(p) >1)
    stop("All elements of p must lie in (0,1]!")
  if(type=="two.sided")
    z <- qnorm(p/2, lower.tail=FALSE)
  if(type=="one.sided")
    z <- qnorm(p, lower.tail=FALSE)
  return(z)
}
