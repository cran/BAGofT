#' @export
#######################################################
# 'testModel' function for binomial regression
#######################################################
testGlmBi <- function(formula, link){
  testModel <- function(Train.data, Validation.data){
    modT <- stats :: glm(formula , family = stats :: binomial(link = link), Train.data)
    # prediction on the training set
    predT <- stats :: predict(modT, newdata = Train.data
                              , type = "response", se.fit = TRUE)$fit
    # prediction on the test set
    predE <- stats :: predict(modT, newdata = Validation.data
                              , type = "response", se.fit = TRUE)$fit

    # calculate the Pearson residual
    res <- stats :: resid(modT, type = "pearson")

    # obtain the response name
    Rsp <- as.character(formula)[2]
    return(list(predT = predT, predE = predE, res = res, Rsp = Rsp))

  }
  return(testModel)
}
