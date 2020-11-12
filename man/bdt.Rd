\name{bdt}
\alias{bdt}
\title{Bootstrap Diagnostic Tool}
\description{bdt is to compute the bias of three ate estimators: inverse probability treatment weighting (IPTW), augmented inverse probability treatment weighting (AIPTW) and targeted maximum likelihood estimator (TMLE) and their corresponding coverage rates. This function returns the propensity scores for all subjects and for subgroups with A=1 and A=0 as well. The \code{bdt} function is generally called with arguments \code{(Y,A,W)}, where \code{Y} is a continuous or binary outcome variable, \code{A} is a binary treatment variable, (\code{A=1} for treatment, \code{A=0} for control), and \code{W} is a matrix or dataframe of baseline covariates. It may contain factors. The number of replications is identified using the \code{M} argument.
}

\usage{bdt(Y, A, W, outcome_type, M, gbound = 0.025, gGLM, Qform = NULL,
gform = NULL, SL.library = NULL)}

\arguments{
  \item{Y}{continuous or binary outcome variable}
  \item{A}{binary treatment indicator, \code{1} - treatment, \code{0} - control}
  \item{W}{vector, matrix, or dataframe containing baseline covariates. It may contain factors.}
  \item{M}{the number of replications to generate bootstrap datasets}
  \item{outcome_type}{"continuous" or "binary"}
  \item{Qform}{regression formula of outcome}
  \item{gGLM}{if TRUE, GLM is used to estimate g; otherwise call \code{SuperLearner}}
  \item{gform}{optional regression formula of g; default value is "A~W"}
  \item{gbound}{value between (0,1) for truncation of predicted probabilities; default value is 0.025}
  \item{SL.library}{vector of prediction algorithms for data adaptive estimation of g}
}

\value{
  \item{true_Effect}{the true ate for bootstrap datasets}
  \item{bias_IPTW}{the IPTW bias ate of bootstrap data based on the true effect}
  \item{bias_AIPTW}{the AIPTW bias ate of bootstrap data based on the true effect}
  \item{bias_TMLE}{the TMLE bias ate by calling \code{tmle} of bootstrap data based on the true effect }
  \item{cov_IPTW}{the coverage rate over \code{M} replications, standard error is computed by calling \code{sandwich}}
  \item{cov_AIPTW}{the coverage rate over \code{M} replications, standard error is computed by influence function}
  \item{cov_TMLE}{the coverage rate over \code{M} replications, standard error is computed by calling \code{tmle}}
  \item{ps_GLM}{the propensity scores estimated by GLM for all subjects over \code{M} replications }
  \item{ps_SL}{the propensity scores estimated by calling \code{SuperLearner} for all subjects over \code{M} replications }
}

\references{
1. Bahamyirou A, Blais L, Forget A, Schnitzer ME. (2019), {Understanding and diagnosing the potential for bias when using machine learning methods with doubly robust causal estimators}. \emph{Statistical methods in medical research}, 28(6), 1637-50.
}

\seealso{

  \code{\link{TMLE_ate}}

  \code{\link{IPTW_ate}}

  \code{\link{AIPTW_ate}}

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
  Qform= "Y~A+W1+W2+I1+I2+P1",
  gform = "A~W2+W1+I1+I2+P1", M = 100)
  summary(result1)

# User-supplied SuperLearner function to estimate g

  result2 <- bdt(Y, A, W, outcome_type ="continuous", gGLM = FALSE,
  Qform= "Y~A+W1+W2+I1+I2+P1",
  SL.library = c("SL.gam", "SL.glm"), M = 100)
  summary(result2)
}

