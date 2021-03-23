\name{bdt}
\alias{bdt}
\title{Bootstrap Diagnostic Tool}
\description{bdt is to compute the bias of two ate estimators: augmented inverse probability treatment weighting (AIPTW) and targeted maximum likelihood estimator (TMLE) and their corresponding coverage rates. This function returns the propensity scores for all subjects and for subgroups with A=1 and A=0 as well. The \code{bdt} function is generally called with arguments \code{(Y,A,W)}, where \code{Y} is a continuous or binary outcome variable, \code{A} is a binary treatment variable, (\code{A=1} for treatment, \code{A=0} for control), and \code{W} is a matrix or dataframe of baseline covariates. It may contain factors. The number of replications is identified using the \code{M} argument.
}

\usage{bdt(Y, A, W, outcome_type, M, seed = NULL, gGLM = NULL, gform = NULL, gbound = 0.025, SL.library = NULL, remove_first_dummy = FALSE, remove_most_frequent_dummy = FALSE)}

\arguments{
  \item{Y}{continuous or binary outcome variable}
  \item{A}{binary treatment indicator, \code{1} - treatment, \code{0} - control}
  \item{W}{vector, matrix, or dataframe containing baseline covariates. It may contain factors.}
  \item{outcome_type}{"continuous" or "binary"}
  \item{M}{the number of replications to generate bootstrap datasets}
  \item{seed}{default value is "NULL". It should be defined if user wants to compare estimators on exactly same bootstrap datasets. The length should equal to M+1}
  \item{gGLM}{if TRUE, GLM is used to estimate g; otherwise call \code{SuperLearner}}
  \item{gform}{optional regression formula of g; default value is "A~W"}
  \item{gbound}{value between (0,1) for truncation of predicted probabilities; default value is 0.025}
  \item{SL.library}{vector of prediction algorithms for data adaptive estimation of g}
  \item{remove_first_dummy}{for categorical covariates, if true remove the first dummy of each covariate such that only n-1 dummies remain. This avoids multicollinearity issues in models(default =\code{FALSE}) }
  \item{remove_most_frequent_dummy}{for categorical covariates, if true remove the most frequently observed category such that only n-1 dummies remain. If there is a tie for most frequent, will remove the first (by alphabetical order) category that is tied for most frequent.(default =\code{FALSE}) }
}

\value{
  \item{true_Effect}{the 'true' ates for bootstrap datasets}
  \item{gbound}{the bounds of ps user defined}
  \item{ps_M}{the propensity scores estimated by GLM of SL for all subjects over \code{M} replications}
  \item{seed}{seeds user defined}
  \item{gGLM}{'GLM' or 'SL' method user defined to estimate ps}
  \item{bias_AIPTW}{the AIPTW bias ate of bootstrap data based on the true effect}
  \item{bias_TMLE}{the TMLE bias ate of bootstrap data based on the true effect }
  \item{cov_AIPTW}{the coverage rate over \code{M} replications, standard error is computed by influence function}
  \item{cov_TMLE}{the coverage rate over \code{M} replications, standard error is computed by influence function}
}

\references{
1. Bahamyirou A, Blais L, Forget A, Schnitzer ME. (2019), {Understanding and diagnosing the potential for bias when using machine learning methods with doubly robust causal estimators}. \emph{Statistical methods in medical research}, 28(6), 1637-50.
}

\seealso{

  \code{\link{summary.bdt}}

  \code{\link{plot.bdt}}
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
  A <- rbinom(n, 1, plogis(0.2+ W[,"W1"]+0.3*W[,"I1"]+W[,"W1"]*W[,"I1"]-0.2*(W[,"W2"]+W[,"I2"])^2 ))
  Y <-  1 + A + W[,"W1"] + 2*W[,"W2"] + 0.5*(W[,"W1"]+W[,"P1"])^2 + rnorm(n)

  result1 <- bdt(Y, A, W, outcome_type = "continuous", gGLM = TRUE,
  gform = "A~W2+W1+I1+I2+P1", M = 60)
  summary(result1)

# User-supplied SuperLearner function to estimate g

  result2 <- bdt(Y, A, W, outcome_type ="continuous", gGLM = FALSE,
  SL.library = c("SL.glm.interaction"), M = 60)
  summary(result2)
}

