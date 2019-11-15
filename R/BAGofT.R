#' @export

# load relevant functions
#source(chi_sel.R)
#source(chi_sel2.R)
#source(splitFUN.R)
#source(groupFUN.R)
#source(splitFUN_P.R)
#source(groupFUN_P.R)
#source(BAGofT_sin.R)

################################################################
#The function for BAGofT test.
################################################################

# nsplits: number of splits
# splitMeth : xqt: use the quantiles of covariates to partition
#            p_fit: use quantiles of fitted probability from the model to partition
#            neru_fit: use quantiles of fitted probability from neural network to partition


BAGofT <- function(formula, data, link = "logit", Ctv = NULL, Dsv = NULL,
                   g = 5, nsplits = 100, spp = 1/2.3, min.Obsr = 10,
                   adj = TRUE, partition.Method = "xqt"){

  ######################################################
  # rename the variables
  fm <- formula
  datset <- data
  lk <- link
  spliDat <- list()
  qtv <- min.Obsr
  splitMeth <- partition.Method
  ######################################################

  for (j in c(1:nsplits)){
    message(paste("split: ", j))
    spliDat[[j]] <- BAGofT_sin(datset = datset, spp = spp, g = g, fm = fm, adj = adj,
                            Ctv = Ctv, Dsv = Dsv, qtv = qtv, lk = lk, splitMeth = splitMeth)
  }
  pvdat <- unlist(lapply(c(1:(nsplits )), function(x) spliDat[[x]]$p.value))
  chidat <- unlist(lapply(c(1:(nsplits)), function(x) spliDat[[x]]$chisq))

  if (splitMeth == "xqt"){
    # in each split, the count of the appearance of each partition covariates
    # used to divide all the sets
    maxgpCtList <- lapply(c(1:(nsplits)), function(x) spliDat[[x]]$maxgpCt)
    # in each split, the count of the appearance of each partition covariates
    # used to divide the set with largest contribution
    allgpCtList <- lapply(c(1:(nsplits)), function(x) spliDat[[x]]$allgpCt)

    # get the aggregated results
    maxgpCtList_Sum <- Reduce(`+`, maxgpCtList)
    allgpCtList_Sum <- Reduce(`+`, allgpCtList)

    nr2 <- nrow(datset)

    # statistic and p values
    BaGofTADstat <- stats :: median(pvdat)
    BaGofTADpval1 <- stats :: pbeta(BaGofTADstat, 50, 50)
    BaGofTADpval2 <- stats :: pnorm(BaGofTADstat, 0.5, sqrt(1/(12 * 100)))
    # using the mean of chisquare statistics
    BaGofTADstat3 <- mean(chidat)
    BaGofTADpval3 <- stats :: pnorm(BaGofTADstat3, 5, sqrt(g*2/nsplits), lower.tail = FALSE)
    message(paste("\n\n\n\n\nMultiple splitting BAGofT\n",
                  "\n Training set size: ",  floor(nr2- g*nr2^spp),
                  "\n Test set size: ",  nr2 - floor(nr2- g*nr2^spp),
                  "\n Number of groups: ", g,
                  "\n Number of splits: ", nsplits,
                  "\n Final sample correction: ", adj,
                  "\n Test statistic value: ", BaGofTADstat,
                  "\n p.value:", BaGofTADpval2))
    return(invisible(list(
      p.value = BaGofTADpval2,
      test.stat = BaGofTADstat,
      p.dat = pvdat,
      chisq.dat = chidat,
      p.value2 = BaGofTADpval1,
      test.stat3 = BaGofTADstat3,
      p.value3 = BaGofTADpval3,
      maxgpCtList_Sum = maxgpCtList_Sum,
      allgpCtList_Sum = allgpCtList_Sum,
      singleSplit.results = spliDat)))
  }else{
    nr2 <- nrow(datset)
    BaGofTADstat <- stats :: median(pvdat)
    BaGofTADpval1 <- stats :: pbeta(BaGofTADstat, 50, 50)
    BaGofTADpval2 <- stats :: pnorm(BaGofTADstat, 0.5, sqrt(1/(12 * 100)))
    # using the mean of chisquare statistics
    BaGofTADstat3 <- mean(chidat)
    BaGofTADpval3 <- stats :: pnorm(BaGofTADstat3, 5, sqrt(g*2/nsplits), lower.tail = FALSE)
    message(paste("\n\n\n\n\nMultiple splitting BAGofT\n",
                  "\n Training set size: ",  floor(nr2- g*nr2^spp),
                  "\n Test set size: ",  nr2 - floor(nr2- g*nr2^spp),
                  "\n Number of groups: ", g,
                  "\n Number of splits: ", nsplits,
                  "\n Final sample correction: ", adj,
                  "\n Test statistic value: ", BaGofTADstat,
                  "\n p.value:", BaGofTADpval2))
    return(invisible(list(
      p.value = BaGofTADpval2,
      test.stat = BaGofTADstat,
      p.dat = pvdat,
      chisq.dat = chidat,
      p.value2 = BaGofTADpval1,
      test.stat3 = BaGofTADstat3,
      p.value3 = BaGofTADpval3,
      singleSplit.results = spliDat)))
  }



}
