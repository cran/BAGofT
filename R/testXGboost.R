#' @export
#######################################################
# 'testModel' function for penalized logistic regression
#######################################################
testXGboost <- function(formula, params = list(), nrounds = 25){
  # return a function of train data and test data
  testModel <- function(Train.data, Validation.data){
    # obtain the response name
    Rsp <- as.character(formula)[2]
    # regressor data
    XmatT <- stats::model.matrix(formula,  Train.data)[,-1]
    XmatE <- stats::model.matrix(formula,  Validation.data)[,-1]
    # fit xgboost
    xgModT <- xgboost :: xgboost(data = XmatT, label = Train.data[, Rsp], params = params, nrounds = nrounds, objective = "binary:logistic", verbose = 0)

    #predict on the test set
    predE <- stats :: predict(xgModT, XmatE)

    #predict on the training set
    predT <- stats :: predict(xgModT, XmatT)

    # calculate the Pearson residual
    res <-   (Train.data[, Rsp] -  predT)/sqrt(predT * (1 - predT ))

    return(list(predT = predT, predE = predE, res = res, Rsp = Rsp))
  }
  return(testModel)
}
