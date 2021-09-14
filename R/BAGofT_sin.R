################################################################
#The function to calculate the chi-square statistic
# value from a single split BAGofT
# previous version fits residual by random forest,
# current version fits Pearson residual
################################################################

BAGofT_sin <- function(testModel, parFun,
                       datset, ne){
  # number of minimum observations in a group
  # nmin <- ceiling(sqrt(ne))

  # number of rows in the dataset
  nr <- nrow(datset)
  # obtain the training set size
  nt <- nr - ne
  # the indices for training set observations
  trainIn <- sample(c(1 : nr), nt)

  #split the data
  datT <- datset[trainIn, ]
  datE <- datset[-trainIn, ]
  # fit the model to test by training data
  testMod <- testModel(Train.data = datT, Validation.data = datE)

  # obtain adaptive partition result from parFun
  par <- parFun(Rsp = testMod$Rsp, predT = testMod$predT, res = testMod$res,
                Train.data = datT, Validation.data = datE)
  #calculate the number of groups left
  ngp <- length(levels(par$gup))

  #########calculate the difference in each group
  dif <- abs(stats :: xtabs(testMod$predE - datE[,testMod$Rsp] ~ par$gup))
  #calculate the denominator in each group
  den <- stats :: xtabs(testMod$predE * (1 - testMod$predE) ~ par$gup)

  #########calculate the chisquare sum
  contri <- (dif)^2/den
  chisq <- sum(contri)

  #calculate test statistic (p value).
  P = 1 - stats :: pchisq(chisq, ngp)
  #pass values to list gls
  gls <- list(chisq = chisq, p.value = P, ngp = ngp, contri = contri, parRes = par$parRes)


  return(gls)
}





########################################################################################
########################################################################################

