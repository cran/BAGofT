################################################################
#Chooses the best partition for a single discrete covariates
#returns the partition that gives the largest "B" value.
################################################################
#x: the covariate to split
#y: the response value
#fp: fitted probability
#mnSt: smallest possible number of observations in
#each group


chi_sel2 <- function(x, y, fp, mnSt) {
  #convert to factor
  x <- as.factor(x)
  #drop the levels with 0 observations
  x <- droplevels(x)
  #get the current column
  vec1 <- levels(x)
  if (length(vec1) <= 1){
    #if there is less than one category
    #assign TRUE to NoDsv
    NoDsv <- TRUE
    result <- list(NoDsv = NoDsv)
  }else{
    #how many selection types: select 1, 2, 3,... for the first part
    #and the n - 1, n - 2, ... for the other part
    cbHalf <- floor(length(vec1)/2)
    combSelList <- list()
    for(cbI in c(1 : cbHalf)){
      combSel <- utils :: combn(seq_along(vec1), cbI)
      
      combSelList <- c(combSelList,
                       lapply(c(1:ncol(combSel)), function(z) combSel[, z]) )
      
    }
    
    #minmum number of observations in each group
    mincount <- unlist(lapply(combSelList, function(z) min(table(x %in% vec1[1]))) )
    
    spSchemeIn <- which(mincount > mnSt)
    if (length(spSchemeIn) == 0){
      #if all of the possible splits do not satisfy the minimum split criterion,
      #assign TRUE to NoDsv
      NoDsv <- TRUE
      result <- list(NoDsv = NoDsv)
    }else{
      #else we may split the discrete covariates
      NoDsv <- FALSE
      ##################################################
      #calculate standardized difference
      ##################################################
      
      #total number of kinds of splitting schemes selected
      spSelNum <- length(spSchemeIn)
      #splitting indices for the first group
      spSel <-   lapply(spSchemeIn, function(z) combSelList [[z]])
      
      #variable to store standardized difference
      chival <- numeric(spSelNum)
      
      for(spSI in c(1 : spSelNum)){
        #index for group 1
        gp1In <- which(x == vec1[spSel[[spSI]]])
        #standardized difference
        chival[spSI] <- (sum(y[gp1In] - fp[gp1In]) )^2/
          sum(fp[gp1In] * (1 - fp[gp1In]) )  +
          (sum(y[-gp1In] - fp[-gp1In]) )^2/
          sum(fp[-gp1In] * (1 - fp[-gp1In]) )
      }
      
      #the point choose to split
      ch <- which.max(chival)
      gp1 <- vec1[spSel[[ch]]]
      gp2 <- vec1[ - spSel[[ch]]]
      result <- list(NoDsv = NoDsv,
                     chival = chival[ch],
                     gp1 = gp1, gp2 = gp2)
    }
  }
  
  
  return(result)
}
