#' @export
#######################################################
# 'testModel' function for penalized logistic regression
#######################################################
testGlmnet <- function(formula, alpha = 1){
  # return a function of train data and test data
  testModel <- function(Train.data, Validation.data){
    # obtain the response name
    Rsp <- as.character(formula)[2]
    # regressor data
    XmatT <- stats::model.matrix(formula,  Train.data)[,-1]
    XmatE <- stats::model.matrix(formula,  Validation.data)[,-1]
    # fit lasso regression
    lasso_cvlamT <- glmnet :: cv.glmnet(XmatT, Train.data[, Rsp], family = "binomial", alpha = alpha)$lambda.min
    lassoModT <- glmnet :: glmnet(XmatT, Train.data[, Rsp], family = "binomial", alpha = alpha, lambda = lasso_cvlamT)

    #predict on the test set
    predE <- stats :: predict(lassoModT, XmatE, type=c("response") )

    #predict on the training set
    predT <- stats :: predict(lassoModT, XmatT, type=c("response") )

    # calculate the Pearson residual
    res <-   (Train.data[, Rsp] -  predT)/sqrt(predT * (1 - predT ))

    return(list(predT = predT, predE = predE, res = res, Rsp = Rsp))
  }
  return(testModel)
}
