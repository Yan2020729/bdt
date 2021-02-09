\name{ps}
\alias{ps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Propensity scores}
\description{
Probabilities of being treated and not being treated for all subjects and for subgroups \code{A=1} and \code{A=0}. Estimates could be computing by using two choices of GLM and/or SuperLearner.
}
\usage{
ps(A, W, gform1 = NULL, gform2 = NULL, SL.library1 = NULL, SL.library2 = NULL, gbound = 0, verbose = TRUE, remove_first_dummy = FALSE, remove_most_frequent_dummy = FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{binary treatment indicator, \code{1} - treatment, \code{0} - control}
  \item{W}{vector, matrix, or dataframe containing baseline covariates}
  \item{gform1}{optional glm regression formula 1 of g}
  \item{gform2}{optional glm regression formula 2 of g}
  \item{SL.library1}{vector of prediction algorithms 1 for data adaptive estimation of g}
  \item{SL.library2}{vector of prediction algorithms 2 for data adaptive estimation of g}
  \item{gbound}{bounds used to truncate g with value between (0,1) for truncation of predicted probabilities; default value is 0}
  \item{verbose}{print the fit summary of GLM or SL if it is TRUE (default=\code{TRUE})}
  \item{remove_first_dummy}{for categorical covariates, if true remove the first dummy of each covariate such that only n-1 dummies remain. This avoids multicollinearity issues in models(default =\code{FALSE}) }
  \item{remove_most_frequent_dummy}{for categorical covariates, if true remove the most frequently observed category such that only n-1 dummies remain. If there is a tie for most frequent, will remove the first (by alphabetical order) category that is tied for most frequent.(default =\code{FALSE}) }
}

\value{
  \item{probabilities}{a dataframe with columns refer to Min., 1st Qu., Median, Mean, 3rd Qu., Max values and rows refer to \eqn{P(A=1|W)} and \eqn{P(A=0|W)} for all subjects, \eqn{P(A=1|W)} for subgroups \code{A=1} and \eqn{P(A=0|W)} subgroups with \code{A=0}  }
  \item{fit_summaries}{summaries of glm regression or SuperLearner models}
}
\examples{
# Example
# Continuous outcome
# User-supplied regression formulas to estimate g

  set.seed(1250)
  n=1000
  sigma <- matrix(c(2, 1, 1, 1), ncol=2)
  X <- matrix(rnorm(n*2), ncol = nrow(sigma)) %*% chol(sigma)
  X <- X + matrix(rep(c(.5, 1),each = n), byrow = FALSE, ncol = 2)
  I1 <- rnorm(n,mean = 1, sd = 2)
  I2 <- rnorm(n,mean = 1, sd = 1.9)
  P1 <- rnorm(n,mean = 1, sd = 1.5)
  W <- data.frame(X, I1, I2, P1)
  colnames(W) <- c("W1", "W2", "I1", "I2",  "P1")
  A <- rbinom(n, 1, plogis(0.2+ W[,"W1"]+0.3*W[,"I1"]+W[,"W1"]*W[,"I1"]-0.2*(W[,"W2"]+W[,"I2"])^2))

  ps1 <- ps(A,W,gform1 = "A~W2+W1+I1+I2",gform2 = "A~W1+I1+I2+W1*W2",
  SL.library1 ="SL.glmnet",gbound=0)

# or

  ps2 <- ps(A,W,gform1 = "A~W2+W1+I1+I2",SL.library1 ="SL.glmnet",
  SL.library2 = "SL.gam",gbound = 0, verbose = FALSE)
}

\references{
1. Bahamyirou A, Blais L, Forget A, Schnitzer ME. (2019), {Understanding and diagnosing the potential for bias when using machine learning methods with doubly robust causal estimators}. \emph{Statistical methods in medical research}, 28(6), 1637-50.
}

