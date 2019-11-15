################################################################
#Chooses the best partition for a single continuous covariates
#returns the partition that gives the largest "B" value.
################################################################
#x: the covariate to split
#y: the response value
#fp: fitted probability
#mnSt: smallest possible number of observations in
#each group

chi_sel <- function(x, y, fp, mnSt) {
  ngroup <- floor(length(x)/mnSt)
  stepS <- 1/ngroup
  qutx <- stats :: quantile(x, seq(stepS, (1 - stepS), stepS  ) )
  chival <- numeric(ngroup - 1)
  for (i in c(1 : (ngroup - 1)) ) {
    #current split point
    sp <- qutx[i]
    #calculate the chisquare value using the current splitting point
    chival[i] <- (sum(y[x < sp] - fp[x < sp]) )^2/
      sum(fp[x < sp] * (1 - fp[x < sp]) )  +
      (sum(y[x >= sp] - fp[x >= sp]) )^2/
      sum(fp[x >= sp] * (1 - fp[x >= sp]) )
    
  }
  #the point choose to split
  ch <- which.max(chival)
  result <- c(chival = chival[ch], split = qutx[ch])
  return(result)
}