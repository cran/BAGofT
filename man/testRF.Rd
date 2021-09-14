\name{testRF}
\alias{testRF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Testing random forests
%%  ~~function to do ... ~~
}
\description{\code{testRF} specifies a random forest as the classifier to test. It returns a function that can be taken as the input of \sQuote{testModel}.
}
\usage{
testRF(formula, ntree = 500, mtry = NULL, maxnodes = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
  an object of class \code{"\link[stats]{formula}"} (or one that
    can be coerced to that class): a symbolic description of the
    model to test.
%%     ~~Describe \code{formula} here~~
}
  \item{ntree}{number of trees to grow. The default is 500.
%%     ~~Describe \code{link} here~~
}
  \item{mtry}{number of variables randomly sampled as candidates at each split. The default value is sqrt(p) where p is the number of covariates.
%%     ~~Describe \code{link} here~~
}
  \item{maxnodes}{maximum number of terminal nodes trees in the forest can have.
%%     ~~Describe \code{link} here~~
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
\dontrun{
###################################################
# Generate a sample dataset.
###################################################
# set the random seed
set.seed(20)
# set the number of observations
n <- 200
# set the number of covariates
p <- 20

# generate covariates data
Xdat <- matrix(runif((n*p), -5,5), nrow = n, ncol = p)
colnames(Xdat) <- paste("x", c(1:p), sep = "")

# generate random coefficients
betaVec <- rnorm(6)
# calculate the linear predictor data
lindat <-  3 * (Xdat[,1] < 2 & Xdat[,1] > -2) + -3 * (Xdat[,1] > 2 | Xdat[,1] < -2) +
  0.5 * (Xdat[,2] + Xdat[, 3] + Xdat[,4] + Xdat[, 5])
# calculate the probabilities
pdat <- 1/(1 + exp(-lindat))

# generate the response data
ydat <- sapply(pdat, function(x) stats :: rbinom(1, 1, x))

# generate the dataset
dat <- data.frame(y = ydat, Xdat)

###################################################
# Obtain the testing result
###################################################

# 50 percent training set
testRes1 <- BAGofT(testModel = testRF(formula = y ~.),
                  data = dat,
                  ne = n*0.5,
                  nsplits = 20,
                  nsim = 40)
# 75 percent training set
testRes2 <- BAGofT(testModel = testRF(formula = y ~.),
                   data = dat,
                   ne = n*0.75,
                   nsplits = 20,
                   nsim = 40)
# 90 percent training set
testRes3 <- BAGofT(testModel = testRF(formula = y ~.),
                   data = dat,
                   ne = n*0.9,
                   nsplits = 20,
                   nsim = 40)

# print the testing result.
print(c(testRes1$p.value, testRes2$p.value, testRes3$p.value))
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
