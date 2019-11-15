################################################################
#Based on the selected partition from splitFUN, group the data
#in the test set datE. Return the group indices for each observation
################################################################


#spR : result from splitFUN
#datE : test data
groupFUN <- function(spR, datE){
  leafs <- spR[which(spR$TERMINAL == "LEAF" | spR$TERMINAL == "SPLIT"), ]
  #if there is one group only,
  if(nrow(leafs) == 1 & is.na(leafs$FILTER[1])){
    group <- rep(1, nrow(datE))
  }else{
    group <- c()
    for (i in seq_len(nrow(leafs))) {
      # extract index
      datE$daind <- c(1 : nrow(datE))
      #get the index of test data observations
      #that fall inside ith group
      #if none of them is in the ith group
      #ind = interger(0) and nothing will be assigned to group
      ind <- subset(datE, eval(parse(
        text = leafs[i, "FILTER"])))$daind
      # estimator is the mean y value of the leaf
      group[ind] <- i
    }
  }
  
  #returns the group index for the test data
  return(as.factor(group))
}