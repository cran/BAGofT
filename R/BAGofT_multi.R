################################################################
#The function to calculate the mean of chi-square statistics
 # from multiple splits.
################################################################

# nsplits: number of splits


BAGofT_multi <- function(testModel, parFun,
                   data,  nsplits, ne){
  spliDat <- list()
  for (j in c(1:nsplits)){
    spliDat[[j]] <- BAGofT_sin(testModel = testModel,
                               parFun = parFun,
                               datset = data,
                               ne = ne)

  }
  # pvalues from multiple splits
  pvdat <- unlist(lapply(c(1:nsplits), function(x) spliDat[[x]]$p.value))

  return( list( meanPv = mean(pvdat),
                medianPv = stats::median(pvdat),
                minPv = min(pvdat),
                spliDat = spliDat )  )
}
