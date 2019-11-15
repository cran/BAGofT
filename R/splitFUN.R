################################################################
#Chooses the best partition based on the continuous covariates
# in ctv and discrete covariates in dsv
################################################################

#rsp: Name of the response variable
# ctv: A character vector that contains the names of
#the continuous variables to partition
# dsv: A character vector that contains the names of
#the categorical variables to partition
# data: Should be the training data to split,
#it contains the response and all covariates needed
# minN: Minimum number of observations in a group
# maxg: Maximum number of groups
# lfit: Fitted probabilites on the training set
#from logistic regression
splitFUN <- function(rsp, ctv, dsv, data, mN, maxg, lfit) {
  #whether these is discrete covariates to partition
  whePdsv <- !is.null(dsv)
  
  #number of observations
  nobs <- nrow(data)
  
  # extract target
  yvec <- data[, rsp]
  
  #logistic fit object
  data$lfit <- lfit
  
  
  #formula for continuous covariates
  formula1 <- stats :: as.formula(paste(rsp, " ~ ", paste(c("-1", ctv),  collapse = " + "), sep = ""))
  
  
  # get the data for continuous covariates. In this way, we allow
  #the continuous covariates to be transformed from the variable
  #in the original data set.
  Xmat1 <- stats :: model.matrix(formula1, data)
  
  
  
  # get the data for continuous covariates. In this way, we allow
  #the continuous covariates to be transformed from the variable
  #in the original data set.
  if(whePdsv){
    Xmat2 <- data[, dsv]
    if(length(dsv == 1)){
      Xmat2 <- data.frame(Xmat2)
      names(Xmat2) <- dsv
    }
  }
  
  
  
  
  
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
        this_data <- subset(data, eval(parse(text = sp_info[j, "FILTER"])))
        # update the design matrix for the subset
        lfit <- this_data$lfit
        this_data$lfit <- NULL
        # get the data for continuous covariates
        Xmat1 <- stats :: model.matrix(formula1, this_data)
        
        
        #get the data for discrete covariates
        if(whePdsv){
          Xmat2 <- this_data[, dsv]
          if(length(dsv == 1)){
            Xmat2 <- data.frame(Xmat2)
            names(Xmat2) <- dsv
          }
        }
        
        # get the data for response
        yvec <- this_data[, rsp]
        
      } else {
        #else, the data is the initial one, no need to update
        this_data <- data
      }
      
      
      
      
      #if none of the covariates can be partitioned, take this
      #node as a leaf and continue to split other node marked
      #as split
      if(nrow(Xmat1) <=  2 *  mN ){
        # overwrite state of current node
        sp_info[j, "TERMINAL"] <- "LEAF"
      }else{
        # overwrite state of current node
        sp_info[j, "TERMINAL"] <- "PARENT"
        
        ###########################split continuous covariates##############
        # estimate splitting criteria to continuous covariates
        splitting1 <- apply(Xmat1,  MARGIN = 2, FUN = chi_sel, y = yvec, fp = lfit,
                            mnSt = mN)
        # get the maximum splitting point
        tmp_splitter <- which.max(splitting1[1,])
        
        
        
        if(whePdsv){
          ###########################split categorical covariates###############
          # estimate splitting criteria to discrete covariates
          splitting2 <- apply(Xmat2,  MARGIN = 2, FUN = chi_sel2, y = yvec, fp = lfit,
                              mnSt = mN)
          
          #find out the covariate that have enough observations to split
          NoDsvVec <- unlist(lapply(c(1 : ncol(Xmat2)), function(z) splitting2[[z]]$NoDsv ) )
          
          #drop the covariates that are not able to split
          covRemain <- which(NoDsvVec == FALSE)
          
          #if all of them are dropped, then we won't partition
          #categorical covariates
          if(length(covRemain) == 0){
            useCatcov <- FALSE
          }else{
            #find out the covariate that give the spliting with largest difference
            dsvChisq <- unlist(lapply(covRemain, function(z) splitting2[[z]]$chival ) )
            
            #best cut by dsv
            dsvSel <- covRemain[which.max(dsvChisq)]
            
            #largest difference from discrete covariate
            disBest <- splitting2[[dsvSel]]$chival
            #largest difference from continuous covariate
            contBest <- splitting1[1, tmp_splitter]
            
            if(disBest > contBest){
              useCatcov <- TRUE
            }else{
              useCatcov <- FALSE
            }
            
          }
          
          
        }
        ################################################################
        
        #if we do not have categorical covariates, then do not use categorical
        #covariates
        if(!whePdsv){
          useCatcov <- FALSE
        }
        #update the  nodes
        if(!useCatcov){
          #if we do not use categorical covariates, then use continuous
          #covariates to split
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
          
          
          
        }else{
          #else, we select the best split in categorical covariates
          
          #group1 levels
          group1String <-  paste("c(",paste(paste("'", splitting2[[dsvSel]]$gp1,"'", sep = ""), collapse = "," ), ")", sep = "")
          
          #group2 levels
          group2String <-  paste("c(",paste(paste("'", splitting2[[dsvSel]]$gp2,"'", sep = ""), collapse = "," ),")",  sep = "")
          
          
          tmp_filter <-   c(paste(names(splitting2)[dsvSel], "%in%",  group1String ),
                            paste(names(splitting2)[dsvSel], "%in%",  group2String ))
          
        }
        
        # define maxnode
        mn <- length(sp_info$NODE)
        
        
        
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
                           FUN = function(i, k) {
                             nrow(subset(x = k, subset = eval(parse(text = i))))
                           },
                           k = this_data)
        
        # create children data frame
        children <- data.frame(NODE = c(mn+1, mn+2),
                               NOBS = tmp_nobs,
                               FILTER = tmp_filter,
                               TERMINAL = rep("SPLIT", 2),
                               row.names = NULL)[split_here,]
        
        
        # bind everything
        sp_info <- rbind(sp_info, children)
        
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
