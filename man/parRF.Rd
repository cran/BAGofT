\name{parRF}
\alias{parRF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Adaptive partition based on random forests
%%  ~~function to do ... ~~
}
\description{\code{parRF} generates an adaptive partition based on the training set data and training set predictions. It controls the group sizes by the covariates data from the validation set.
}
\usage{
parRF(parVar = ".", Kmax = NULL, nmin = NULL, ntree = 60, mtry = NULL, maxnodes = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parVar}{ a character vector that contains the names of the  covariates to generate the adaptive partition. The default is taking all the variables except the response from the \sQuote{data}.
%%     ~~Describe \code{formula} here~~
}
  \item{Kmax}{ the maximum number of groups. The default is floor(nrow(Test.data)/nmin).
%%     ~~Describe \code{formula} here~~
}
  \item{nmin}{ a numerical vector of the training set Pearson residuals from the classifier to test. The default is ceiling(sqrt(nrow(Validation.data))).
%%     ~~Describe \code{formula} here~~
}
  \item{ntree}{number of trees to grow. The default is 60.
%%     ~~Describe \code{link} here~~
}
  \item{mtry}{number of variables randomly sampled as candidates at each split. The default value is floor( length(parVarNew)/3) where parVarNew is the number of covariates after the preselection.
%%     ~~Describe \code{link} here~~
}

  \item{maxnodes}{maximum number of terminal nodes trees in the forest can have.
%%     ~~Describe \code{link} here~~
}

}
\value{
    \item{gup}{ a factor that contains the grouping result of the validation set data.
  }
    \item{parRes}{ a list that contains the variable importance from the random forest.
  }
}

\references{Zhang, Ding and Yang (2021) "Is a Classification Procedure Good Enough?-A Goodness-of-Fit Assessment Tool for Classification Learning" arXiv preprint 	arXiv:1911.03063v2 (2021).
%%% ~put references to the literature/web site here ~
}
%\author{
%%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
###################################################
# Generate a sample dataset.
###################################################
# set the random seed
set.seed(20)
# set the number of observations
n <- 200

# generate covariates data
x1dat <- runif(n, -3, 3)
x2dat <- rnorm(n, 0, 1)
x3dat <- rchisq(n, 4)

# set coefficients
beta1 <- 1
beta2 <- 1
beta3 <- 1

# calculate the linear predictor data
lindat <- x1dat * beta1 + x2dat * beta2 + x3dat * beta3
# calculate the probabilities by inverse logit link
pdat <- 1/(1 + exp(-lindat))

# generate the response data
ydat <- sapply(pdat, function(x) stats :: rbinom(1, 1, x))

# generate the dataset
dat <- data.frame(y = ydat, x1 = x1dat, x2 = x2dat,
                  x3 = x3dat)

###################################################
# Apply parRF to generate an adaptive partition
###################################################
# number of rows in the dataset
nr <- nrow(dat)
# size of the validation set
ne <- floor(5*nrow(dat)^(1/2))
# obtain the training set size
nt <- nr - ne
# the indices for training set observations
trainIn <- sample(c(1 : nr), nt)

#split the data
datT <- dat[trainIn, ]
datE <- dat[-trainIn, ]
# fit a logistic regression model to test by training data
testModel <- testGlmBi(formula = y ~ x1 + x2 , link = "logit")
# output training set predictions and pearson residuals
testMod <- testModel(Train.data = datT, Validation.data = datE)

# obtain adaptive partition result from parFun
parFun <- parRF(parVar = c("x1", "x2", "x3"))
par <- parFun(Rsp = testMod$Rsp, predT = testMod$predT, res = testMod$res,
              Train.data = datT, Validation.data = datE)

# print the grouping result of the validataion set data
print(par$gup)

# print variable importance from the random forest
print(par$parRes)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
