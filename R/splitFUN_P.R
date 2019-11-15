################################################################
#Chooses the best partition based on the fitted probability
# from one of "xqt", "p_fit", "neu_fit" and "rf_fit".
################################################################

#rsp: Name of the response variable
# phat
# data: Should be the training data to split,
#it contains the response and all covariates needed
# minN: Minimum number of observations in a group
# maxg: Maximum number of groups
# lfit: Fitted probabilites on the training set
#from logistic regression


splitFUN_P <- function(rsp, phat, data, mN, maxg, lfit) {
  nobs <- nrow(data)
  # get the data for continuous covariates. In this way, we allow
  #the continuous covariates to be transformed from the variable
  #in the original data set.
  temp_data <- data.frame(phat = phat)
  
  # extract target
  yvec <- data[, rsp]
  
  # add to the temporary data frame
  temp_data[[rsp]] <- yvec
  
  #formula for continuous covariates
  formula1 <- stats :: as.formula(paste(rsp, " ~ -1+ phat", sep = ""))
  
  Xmat1 <- stats :: model.matrix(formula1, temp_data)
  
  
  
  temp_data$lfit <- lfit
  # initialize while loop
  do_splits <- TRUE
  
  # create output data.frame with splitting rules and observations
  sp_info <- data.frame(NODE = 1, NOBS = nrow(data), FILTER = NA,
                        TERMINAL = "SPLIT",
                        stringsAsFactors = FALSE)
  
  while(do_splits) {
    
    # which parents have to be splitted
    to_calculate <- which(sp_info$TERMINAL == "SPLIT")
    for (j in to_calculate) {
      # handle root node
      if (!is.na(sp_info[j, "FILTER"])) {
        #the data is not the initial data
        # subset data according to the filter
        this_data <- subset(temp_data, eval(parse(text = sp_info[j, "FILTER"])))
        # update the design matrix for the subset
        lfit <- this_data$lfit
        this_data$lfit <- NULL
        # get the data for continuous covariates
        Xmat1 <- stats :: model.matrix(formula1, this_data)
        
        # get the data for response
        yvec <- this_data[, rsp]
        
      } else {
        #else, the data is the initial one, no need to update
        this_data <- temp_data
      }
      
      
      
      
      #if none of the covariates can be partitioned, take this
      #node as a leaf and continue to split other node marked
      #as split
      # if this group doesn't
      # have enough unique values to find all the required quantiles
      # stop splitting
      if(nrow(Xmat1) <=  2 *  mN | length(unique(Xmat1[,1])) < floor(length(Xmat1[,1])/mN)){
        # overwrite state of current node
        sp_info[j, "TERMINAL"] <- "LEAF"
      }else{
        # overwrite state of current node
        sp_info[j, "TERMINAL"] <- "PARENT"
        # estimate splitting criteria to continuous covariates
        splitting1 <- apply(Xmat1,  MARGIN = 2, FUN = chi_sel, y = yvec, fp = lfit,
                            mnSt = mN)
        # sometimes, not enough different phat values causes problem for
        # the splitting, and returns splitting1 with numeric(0)
        # if that happens, make that terminal to be a leaf.
        if (length(splitting1) == 0){
          # overwrite state of current node
          sp_info[j, "TERMINAL"] <- "LEAF"
        }else{
          # get the maximum splitting point
          tmp_splitter <- which.max(splitting1[1,])
          # define maxnode
          mn <- length(sp_info$NODE)
          # paste filter rules
          #if there is only one covariate to be considered,
          #tmp_splitter will be a vector, thus has no name
          #so assign the column name of splitting instead
          if (ncol(splitting1) == 1 ){
            tmp_filter <- c(paste(colnames(splitting1), ">=",
                                  splitting1[2,tmp_splitter]),
                            paste(colnames(splitting1), "<",
                                  splitting1[2,tmp_splitter]))
          }else{
            tmp_filter <- c(paste(names(tmp_splitter), ">=",
                                  splitting1[2,tmp_splitter]),
                            paste(names(tmp_splitter), "<",
                                  splitting1[2,tmp_splitter]))
          }
          
          
          # Error handling! check if the splitting rule has already been invoked
          #the x in x=y is the argument name, which is different from the
          #previous x.
          # split_here  <- !sapply(tmp_filter,
          #                       FUN = function(x,y) any(grepl(x, x = y)),
          #                      y = sp_info$FILTER)
          
          split_here <- c(TRUE, TRUE)
          # append the splitting rules
          #if previous one is not the initial node, append the tree_info
          #to it.
          if (!is.na(sp_info[j, "FILTER"])) {
            tmp_filter  <- paste(sp_info[j, "FILTER"],
                                 tmp_filter, sep = " & ")
          }
          
          # get the number of observations in current node
          tmp_nobs <- sapply(tmp_filter,
                             FUN = function(i, x) {
                               nrow(subset(x = x, subset = eval(parse(text = i))))
                             },
                             x = this_data)
          
          # create children data frame
          children <- data.frame(NODE = c(mn+1, mn+2),
                                 NOBS = tmp_nobs,
                                 FILTER = tmp_filter,
                                 TERMINAL = rep("SPLIT", 2),
                                 row.names = NULL)[split_here,]
          
          
          # bind everything
          sp_info <- rbind(sp_info, children)
          
        }
        
      }
      
      
      #number of groups equals the product of the
      #number of splits in each value plus one
      ng <- sum(sp_info$TERMINAL == "LEAF" | sp_info$TERMINAL == "SPLIT")
      # check if there are any open splits left
      
      do_splits <- any(sp_info$TERMINAL == "SPLIT") & (ng < maxg)
      if (!do_splits) {
        break
      }
    }
  }
  #return the information about which point to split
  return(sp_info)
}
