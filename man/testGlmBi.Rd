\name{testGlmBi}
\alias{testGlmBi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Testing binomial regressions
%%  ~~function to do ... ~~
}
\description{\code{testGlmBi} specifies a binomial regression as the classifier to test. It returns a function that can be taken as the input of \sQuote{testModel}.
}
\usage{
testGlmBi(formula, link)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
  an object of class \code{"\link[stats]{formula}"} (or one that
    can be coerced to that class): a symbolic description of the
    model to test.
%%     ~~Describe \code{formula} here~~
}
  \item{link}{a specification for the model link function. Can be one of "logit", "probit", "cloglog".
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
# Obtain the testing result
###################################################
# Test a logistic regression that misses 'x3'. The partition
# variables are 'x1', 'x2', and 'x3'.
testRes <- BAGofT(testModel =testGlmBi(formula = y ~ x1 + x2 , link = "logit"),
       parFun = parRF(parVar = c("x1", "x2", "x3")),
       data = dat)

# the bootstrap p-value is 0. Therefore, the test is rejected
print(testRes$p.value)

# the variable importance from the adaptive partion shows that x3 is likely
# to be the reason for the overfitting (,which is correct since the formula
# fm misses the x3).
print(VarImp(testRes))
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
