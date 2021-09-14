#' @export
# calculates variable importance from test output

VarImp <- function(TestRes){
  viList <- TestRes$singleSplit.results
  Var.imp <- apply(simplify2array(lapply(c(1:length( viList)), function(x){viList[[x]]$parRes$Var.imp})), c(1,2), mean)

  if ("preVar.imp" %in% names(TestRes$singleSplit.results[[1]]$parRes) ) {
    preVar.imp <- apply(simplify2array(lapply(c(1:length( viList)), function(x){viList[[x]]$parRes$preVar.imp})), c(1,2), mean)
    return(list(Var.imp = Var.imp, preVar.imp = preVar.imp))
    }else{
    return(list(Var.imp = Var.imp))
  }

  }



