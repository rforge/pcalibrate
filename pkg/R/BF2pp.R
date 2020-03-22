
# transformation of BF to posterior probability (of the hypothesis in the numerator of the bf)
BF2pp <- function(BF, prior.prob=0.5){
  prior.odds <- prior.prob/(1-prior.prob)
  log.post.odds <- log(prior.odds) + log(BF)
  post.prob <- exp(log.post.odds)/(1+ exp(log.post.odds))
  return (post.prob)
}