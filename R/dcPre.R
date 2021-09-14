# preselection by distance correlation

dcPre <- function( npreSel = 5,
                  type = "V"){
  
  return(function(datRf, parVar){
    # count the number of variables to partition
    nParV <- if(!identical(parVar, ".")){length(parVar)}else{ncol(datRf) - 1}

    if (nParV > npreSel){
      
      preSelected <- TRUE
      
      datRf_temp <- datRf
      datRf_temp$res <- NULL
      # calculate the distance correlation
      dcRes <- t(dcov :: mdcor(datRf$res, datRf_temp))
      rownames(dcRes) <- colnames(datRf_temp)
      # select Kmax partition variables with the largest variable importance
      parVarNew <- colnames(datRf_temp)[order(-as.numeric(dcRes))[c(1: npreSel)]]
      return(list(preSelected = preSelected, parVarNew = parVarNew, VI = dcRes))
      
    }else{
      preSelected = FALSE
      parVarNew <- parVar
      return(list(preSelected = preSelected, parVarNew = parVarNew))
      
    }
    
  })
}