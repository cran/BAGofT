\name{BAGofT}
\alias{BAGofT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A Binary Regression Adaptive Goodness-of-fit Test (BAGofT)
%%  ~~function to do ... ~~
}
\description{\code{BAGofT} is used to test the goodness-of-fit of binary classifiers.
The test statistic is constructed based on the results from multiple splittings.
In each split, the test first
splits the data into a training set and a validation set. Then,
it adaptively obtains a partition based on the training set and performs a goodness-of-fit test on the validation set. Details can be found in Zhang, Ding and Yang (2021).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BAGofT(testModel, parFun = parRF(), data, nsplits = 100,
ne = floor(5*nrow(data)^(1/2)), nsim = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{testModel}{ a function that generates predicted results from the classifier to test. Details can be found in \code{"\link[BAGofT]{testGlmBi}"} for binomial regression, \code{"\link[BAGofT]{testGlmnet}"} for penalized logistic regression, \code{"\link[BAGofT]{testRF}"} for random forest, and \code{"\link[BAGofT]{testXGboost}"} for XGboost.
%%     ~~Describe \code{formula} here~~
}
  \item{parFun}{ a function that generates the adaptive partition. The default is \sQuote{parRF()} that generates a partition by random forest. More information can be found in \code{"\link[BAGofT]{parRF}"}.
%%     ~~Describe \code{formula} here~~
}
 \item{data}{a data frame containing  the response and covariates used in the model together with the other covariates not in the model but considered used to generate the partition.
%%     ~~Describe \code{data} here~~
}
  \item{nsplits}{number of splits. The default is 100.
%%     ~~Describe \code{nsplits} here~~
}
  \item{ne}{the size of the validation set. The default is floor(5*nrow(data)^(1/2)).
%%     ~~Describe \code{nsplits} here~~
}
  \item{nsim}{the number of simulated datasets to calculate the bootstrap \eqn{p}-value.
%%     ~~Describe \code{nsplits} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above %~~
%}
\value{
    \item{p.value}{the bootstrap \eqn{p}-value of the BAGofT test statistic (which combines the results from multiple splitting by taking the average).
  }
    \item{p.value2}{the bootstrap \eqn{p}-value from an alternative version of the BAGofT test statistic (which combines the results from multiple splitting by taking the sample median).
  }
    \item{p.value3}{the bootstrap \eqn{p}-value from an alternative version of the BAGofT test statistic (which combines the results from multiple splitting by taking the minimum).
  }
    \item{pmean}{the BAGofT test statistic (which combines the results from multiple splitting by taking the average).
  }
    \item{pmedian}{an alternative BAGofT test statistic (which combines the results from multiple splitting by taking the sample median).
  }
    \item{pmin}{an alternative BAGofT test statistic (which combines the results from multiple splitting by taking the minimum).
  }
    \item{simRes}{a list that contains the simulated test statitics used to generate the bootstrap \eqn{p}-values. \sQuote{simRes$pmeanSim}, \sQuote{simRes$pmediansim}, \sQuote{simRes$pmeanSim} corresepond to the three kinds of BAGofT statistics, respectively.
  }
    \item{singleSplit.results}{a list that contains the results from each splitting. Its elements are as follows.

\sQuote{singleSplit.results[[k]]$chisq}: The chi-squared statistic of the BAGofT test from the \eqn{k}th splitting.

\sQuote{singleSplit.results[[k]]$p.value}: The \eqn{p}-value calculated from the chi-squared statistic.

\sQuote{singleSplit.results[[k]]$ngp}: The number of groups chosen by the adaptive partition.

\sQuote{singleSplit.results[[k]]$contri}: The weighted sum of squares from each group.

\sQuote{singleSplit.results[[k]]$parRes}: Variable importance (or other results from custom partition functions) from the adaptive partition.
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

# the variable importance from the adaptive partition shows that x3 is likely
# to be the reason for the overfitting (,which is correct since the formula
# fm misses the x3).
print(VarImp(testRes))
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
