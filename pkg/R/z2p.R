
z2p <- function(z, type="two.sided"){
  if(type=="two.sided")
    p <- 2*pnorm(abs(z), lower.tail=FALSE)
  if(type=="one.sided")
    p <- pnorm(z, lower.tail=FALSE)
  return(p)
}
