#' @export
parRF <- function(parVar = ".",
                  Kmax = NULL,
                  nmin = NULL,
                  ntree = 60,
                  mtry = NULL,
                  maxnodes = NULL
                  ){
  # return a function of response, training set fitted values, training set
  # residual, training data and test data.
  return(function(Rsp, predT, res,
                   Train.data, Validation.data){
    # apply a build-in function to preselect partition covaraites
    # by distance correlation
    preFun = dcPre()
    if (is.null(nmin)){
      nmin <- ceiling(sqrt(nrow(Validation.data)))
    }
    # assign default vallues
    if (is.null(Kmax)){
      Kmax <- floor(nrow(Validation.data)/nmin)
    }

    if (is.null(maxnodes)){
      maxnodes <- if(!identical(parVar, ".")){min(nrow(Train.data), ceiling(5 * length(parVar)) )}else{min(nrow(Train.data), ceiling(5 * ( ncol(Train.data) - 1)) )}
    }

    # training data with phat - y
    datRf <- Train.data[, -which(names(Train.data) == Rsp)]

    # calculate the Pearson residual
    datRf$res <- res

    # preselection result
    preSelRes <- preFun(datRf = datRf, parVar = parVar)
    parVarNew <- preSelRes$parVarNew
    preSelected <- preSelRes$preSelected
    if (is.null(mtry)){
      mtry <- floor( length(parVarNew)/3)
    }
    # partition by random forest based on preselected variables
    formula_rf <- stats :: as.formula(paste("res", " ~ ", paste(parVarNew,  collapse = " + "), sep = ""))

    resRf <- randomForest :: randomForest(formula_rf, data = datRf, ntree = ntree,  maxnodes = maxnodes, mtry = mtry, importance=TRUE)
    # obtain random forest residual prediction on the training set
    trainsetPred <-  stats :: predict(resRf, newdata = Train.data)

    # sometimes prediction on the training set will result in all 1 or 0
    # and grouping by the quantiles of the random forest will be problematic
    # since there is only one unique value.
    # if it happens, take all the observations into one group.
    if (length( unique(round(trainsetPred, digits = 10)) ) == 1){
      #the number of groups left
      ngp <- 1

      gup <- as.factor(rep(1, nrow(Validation.data)) )
    }else{
      # else, divide the groups by the quantiles.

      # adjusted maximum number of groups, in case there are too few unique fitted values
      Kmax_adj <- min(Kmax, length( unique(round(trainsetPred, digits = 10)) ) )
      # from 1 to Kmax, calculate the chi-square value on the training set
      chitrainVec <- numeric(Kmax_adj)
      for (gt in c(1:Kmax_adj)){
        # divide the fitted values on the training set by training set quantiles
        gupt <- cut(round(trainsetPred, digits = 10) ,
                    breaks = stats::quantile(unique(round(trainsetPred, digits = 10)) , probs = seq(0,  1, 1/gt)), include.lowest = TRUE)
        #########calculate the difference in each group on the training set
        dift <- abs(stats :: xtabs(  predT - Train.data[,Rsp] ~ gupt))
        #calculate the denominator in each group
        dent <- stats :: xtabs(predT * (1 - predT) ~ gupt)
        #########calculate the test statistic
        contrit <- (dift)^2/dent
        chitrainVec[gt] <- sum(contrit)
        # if one group in the test set has no observation, the contribution is NaN
        # if that happens fill the value by the value from k-1. The previous
        # k always be selected with high priority than this group
        if(is.nan(sum(contrit))){
          chitrainVec[gt] <- chitrainVec[(gt-1)]
        }
      }

      # Select the number of groups
      chitrainVec2 <- chitrainVec[-1]
      chitrainVec3 <- chitrainVec[-length(chitrainVec)]
      chitrainVec4 <- chitrainVec2 - chitrainVec3
      Ksel <- which.max(chitrainVec4) + 1

      # obtain random forest prediction on the test set
      testsetPred <-  stats :: predict(resRf, newdata = Validation.data)

      # divide the prediction on the test set by training set quantiles
      gup <- cut(testsetPred ,
                 breaks = c(Inf, -Inf, stats::quantile(unique(round(trainsetPred, digits = 10)) , probs = seq(0,  1, 1/Ksel)) ) , include.lowest = TRUE)

      #drop the levels with 0 observations
      gup <- droplevels(gup)
      # if the number of observations in the smallest group is less than
      # nmin, combine it with the second smallest group. Repeat this until
      # all of the groups have size not less than nmin
      freqTab <- table(gup)
      while(min(freqTab) < nmin){
        # combine the smallest level with the second smallest level
        levels(gup)[which(levels(gup) == names(sort(freqTab))[1])] <- names(sort(freqTab))[2]
        # update the frequency table
        freqTab <- table(gup)
      }
    }

    # store partition result
    if (preSelected){
      # store preselection importance
      Pre_importance <- preSelRes$VI
      # create a zero matrix
      selImp <- matrix(0, nrow = nrow(Pre_importance), ncol = 2)
      rownames(selImp) <- rownames(Pre_importance)
      colnames(selImp) <- c("%IncMSE", "IncNodePurity")
      # store the variable importance of selected variables
      matchInd <- sapply(rownames(resRf$importance), function(x) which(rownames(Pre_importance) == x) )
      selImp[matchInd, ] <- resRf$importance

      parRes <- list(Var.imp = selImp,
                     preVar.imp = Pre_importance)
    }else{
      parRes <- list(Var.imp = resRf$importance)
    }

    return(list(gup = gup, parRes = parRes))
  })

}

