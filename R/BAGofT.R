#' @export

################################################################
#The function to calculate the emperical p-value of multiple-split
 # results from simulated dataset
 # If the number of partition variables is larger than Kmax,
 # apply a random forest to preselect Kmax partition variables
################################################################
# g: number of groups from the partition
# nsplits: number of splits applied in the multiple splitting
# ne: test set size
# nsim: number of simulated data sets to calculate emperical p-value
BAGofT <- function(testModel, parFun = parRF(),
                   data, nsplits = 100, ne = floor(5*nrow(data)^(1/2)), nsim = 100){
  testRes <- BAGofT_multi(testModel = testModel, parFun = parFun,
                           data = data,  nsplits = nsplits, ne = ne)



  if (nsim >= 1){
    # simulate data
    pmeansimVec <- numeric(nsim)
    pmediansimVec <- numeric(nsim)
    pminsimVec <- numeric(nsim)

    message("Generating simulated data for empirical p-value")

    # fit the model to test by training data
    # we do no need the output for the test set prediction
    # take a dataset with a single row with all 0s as input
    dataFittemp <-  data
    dataFittemp[1, ]<- 0
    # fit model on the dataset
    modFit <- testModel(Train.data = data, Validation.data = dataFittemp)
    # probability calculated from fitted coefficients
    pdat2 <- modFit$predT
    for (i in c(1:nsim)){
      # process bar
      message(paste("Calculating results from ", i, "th simulated dataset"))

      # randomly generated data from the fitted probabilities
      ydat2 <- sapply(pdat2, function(x) stats::rbinom(1, 1, x))
      dat2 <- data
      Rsp <- modFit$Rsp
      dat2[,Rsp] <- ydat2
      # mean statistic calculated from the simulated data.
      testRes_sim <-  BAGofT_multi(testModel = testModel, parFun = parFun,
                                   data = dat2,  nsplits = nsplits, ne = ne)

      pmeansimVec[i] <- testRes_sim$meanPv
      pmediansimVec[i] <- testRes_sim$medianPv
      pminsimVec[i] <- testRes_sim$minPv
    }
    # calculate empirical p-value from simulated data
    pvalue <- mean(testRes$meanPv > pmeansimVec)
    pvalue2 <- mean(testRes$medianPv > pmediansimVec)
    pvalue3 <- mean(testRes$minPv > pminsimVec)

    message(paste("p-value: ",  pvalue,
                  "Averaged statistic value: ",  testRes$meanPv))
    return(invisible( list(p.value = pvalue,
                           p.value2 = pvalue2,
                           p.value3 = pvalue3,
                           pmean = testRes$meanPv,
                           pmedian = testRes$medianPv,
                           pmin = testRes$minPv,
                           simRes = list(pmeanSim = pmeansimVec,
                                         pmediansim = pmediansimVec,
                                         pminsim = pminsimVec),
                           singleSplit.results = testRes$spliDat )      ))
  }else{
    message(paste("Averaged statistic value: ",  testRes$meanPv))
    return(     invisible( list(
                           pmean = testRes$meanPv,
                           pmedian = testRes$medianPv,
                           pmin = testRes$minPv,
                           singleSplit.results = testRes$spliDat )      ))

  }




}





########################################################################################
########################################################################################

