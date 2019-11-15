################################################################
#The function for single BAGofT test.
################################################################


# ctv: A character vector that contains the names of
#the continuous variables to partition
# dsv: A character vector that contains the names of
#the categorical variables to partition

BAGofT_sin <- function(datset, spp = 1/2.3, g,
                    fm, adj = TRUE, Ctv = NULL, Dsv = NULL,qtv,
                    lk = "logit", splitMeth = "xqt"){
  
  if(is.null(Ctv) & is.null(Dsv) & splitMeth == "xqt"){
    stop("Need to input at least one of ctv and dsv")
  }
  
  nr <- nrow(datset)
  nt <- nr - g*nr^spp
  #minimum number of obserations in a group
  minN <- floor(nr/qtv)
  trainIn <- sample(c(1 : nr), nt)
  #split the data
  datT <- datset[trainIn, ]
  datE <- datset[-trainIn, ]
  #fit the logistic regression on the training set
  modT <- stats :: glm(fm , family = stats :: binomial(link = lk), data=datT)
  
  #predict on the test set
  predM <- stats :: predict(modT, newdata = datE
                            , type = "response", se.fit = TRUE)
  #fitted value
  predE <- predM$fit
  
  
  ###Use splitFUN to find a maximal split##########################
  #get the name of response variable
  Rsp <- as.character(fm)[2]
  
  if (splitMeth == "xqt"){
    spResult <- splitFUN(rsp = Rsp, ctv = Ctv, dsv = Dsv,
                         data = datT, mN = minN,
                         maxg = g, lfit = stats :: fitted(modT))
    # get the grouping results for the test data
    gup <- groupFUN(spR = spResult, datE = datE)
    
    
  }else if (splitMeth == "p_fit"){
    spResult <- splitFUN_P(rsp = Rsp, phat = stats :: fitted(modT),
                           data = datT, mN = minN,
                           maxg = g, lfit = stats :: fitted(modT))
    # get the grouping results for the test data
    gup <- groupFUN_P(spR = spResult, datE = datE, phatE = predE)
  }else if (splitMeth == "neu_fit"){
    # if we want to use the p-value from neural network we need
    # the package nnet.
    #library(nnet)
    
    formula_neu <- stats :: as.formula(paste("factor(", Rsp,")", " ~ ", paste(c(Ctv, Dsv),  collapse = " + "), sep = ""))
    neum <- nnet :: nnet(formula_neu, data=datT, size=5, decay=0.01, maxit=500, trace = FALSE)
    
    
    #predicted value using neuro network
    # training set fitted value
    npreT <- stats :: predict(neum, datT, type="raw")
    # test set predicted value
    npreE <- stats :: predict(neum, datE, type="raw")
    
    spResult <- splitFUN_P(rsp = Rsp, phat = npreT,
                           data = datT, mN = minN,
                           maxg = g, lfit = stats :: fitted(modT))
    # get the grouping results for the test data
    gup <- groupFUN_P(spR = spResult, datE = datE, phatE = npreE)
  }else if (splitMeth == "rf_fit"){
    # if we want to use the p-value from neural network we need
    # the package randomForest.
    #library(randomForest)
    
    formula_rf <- stats :: as.formula(paste("factor(", Rsp,")", " ~ ", paste(c(Ctv, Dsv),  collapse = " + "), sep = ""))
    
    rfm <- randomForest :: randomForest(formula_rf, data = datT, importance=TRUE)
    
    #predicted value using neuro network
    # training set fitted value
    rfpT <- stats :: predict(rfm, datT, type = "prob",
                             norm.votes = TRUE, predict.all = FALSE,
                             proximity = FALSE, nodes = FALSE)[,2]
    # test set predicted value
    rfpE <- stats :: predict(rfm, datE, type = "prob",
                             norm.votes = TRUE, predict.all = FALSE,
                             proximity = FALSE, nodes = FALSE)[,2]
    
    spResult <- splitFUN_P(rsp = Rsp, phat = rfpT ,
                           data = datT, mN = minN,
                           maxg = g, lfit = stats :: fitted(modT))
    # get the grouping results for the test data
    gup <- groupFUN_P(spR = spResult, datE = datE, phatE =  rfpE)
    
  }
  
  #calculate the number of groups left
  ngp <- length(levels(gup))
  
  #########calculate the difference in each group
  dif <- abs(stats :: xtabs(predE - datE[,Rsp] ~ gup))
  #calculate the denominator in each group
  den <- stats :: xtabs(predE * (1 - predE) ~ gup)
  
  #########calculate the test statistic
  contri <- (dif)^2/den
  chisq <- sum(contri)
  
  
  if(adj == TRUE){
    
    ##########get the model object#######################
    sT <- summary(modT)
    #get the data for the current covariate
    datmat <- stats :: model.matrix(fm, datE)
    nc <- ncol(datmat)
    betaDeri <- numeric(nc)
    for(p in c(1:nc)){
      xvec <- datmat[, p]
      
      predL <- stats :: predict(modT, newdata = datE
                                , type = "link")
      
      part1 <- stats :: xtabs(predE*(1 -predE) ~ gup)
      
      part2 <- stats :: xtabs( datE[, Rsp]-predE ~ gup)
      
      #the partial derivative of phat
      if (lk == "logit"){
        deri <- exp(predL)/(1 + exp(predL))^2 * xvec
      }else if(lk == "probit"){
        deri <- stats :: dnorm(predL) * xvec
      }else if(lk == "cloglog"){
        deri <- exp(predL - exp(predL)) * xvec
      }else{
        stop("invalid link")
      }
      
      part3 <-  stats :: xtabs( deri ~ gup)
      
      part4 <- stats :: xtabs( (1 - 2*predE)*deri ~ gup )
      #derivative with respect to the pth coefficient
      betaDeri[p] <- sum((-2 * part2 * part1 * part3  - part4 * part2^2)/part1^2)
      
    }
    
    #calculate estimated standard error
    seHat <- sqrt(as.numeric(t(betaDeri) %*%sT$cov.unscaled %*%   betaDeri))
    
    #calculate the adjusted statistic
    chisq <- max(chisq - seHat * stats :: qnorm(0.95), 0 )
  }
  
  #calculate p value.
  P = 1 - stats :: pchisq(chisq, ngp)
  
  # extract leafs and contribution information
  leafs <- spResult[which(spResult$TERMINAL == "LEAF" | spResult$TERMINAL == "SPLIT"), ]
  
  #library(stringr)
  # which set gives the maximum contribution
  maxgup <- which.max(contri)
  
  # get a vector of variable names to count from the leafs
  if (is.null(Dsv)){
    varlist <- Ctv
  }else{
    varlist <- c(Ctv, Dsv)
  }
  varlist2 <- paste(varlist, " ", sep = "")
  
  # if the split is based on the quantile of covariates,
  # out put the frequency of each covariates used
  if (splitMeth == "xqt"){
    # count the number of appearance of each splitting variable
    # in the group that gives maximum contribution
    maxgpCt <- sapply(varlist2, function(x)stringr :: str_count(leafs$FILTER[maxgup], pattern = x))
    # count the the number of appearance of each splitting variable
    # in all of the groups
    allgpCt <- sapply(varlist2, function(x)stringr :: str_count(paste(leafs$FILTER, collapse = ''), pattern = x))
    
    maxleaf <- leafs$FILTER[maxgup]
    #pass values to list gls
    gls <- list(chisq = chisq, p.value = P, ngp = ngp, leafs = leafs, contri = contri,
                maxgup = maxgup, maxgpCt = maxgpCt, allgpCt = allgpCt,
                maxleaf = maxleaf)
  }else{
    #pass values to list gls
    gls <- list(chisq = chisq, p.value = P, ngp = ngp, leafs = leafs, contri = contri)
  }
  
  return(gls)
}





########################################################################################
########################################################################################

