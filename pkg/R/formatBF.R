
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