\name{summary.bdt}
\alias{summary.bdt}
\alias{print.summary.bdt}
\title{Summarization of the results of a call to the bdt function}
\description{These functions are \link{methods} for class \code{bdt}, \code{summary.bdt} objects}
\usage{
  \method{summary}{bdt}(object, \dots)
  \method{print}{summary.bdt}(x, \dots)
}

\arguments{
  \item{object}{an object of class \code{bdt} for summary functions.}
  \item{x}{an object of class \code{summary.bdt} for print functions.}
  \item{\dots}{currently ignored.}
}

\value{
  \item{true_Effect}{the true ate with respect to bootstrap datasets}
  \item{bias.iptw}{the Min., 1st Qu., Median, Mean, 3rd Qu., Max value of IPTW bias ate over \code{M} replications }
  \item{bias.aiptw}{the Min., 1st Qu., Median, Mean, 3rd Qu., Max value of AIPTW bias ate over \code{M} replications }
  \item{bias.tmle}{the Min., 1st Qu., Median, Mean, 3rd Qu., Max value of TMLE bias ate over \code{M} replications }
  \item{cov.iptw}{the coverage rate over \code{M} replications}
  \item{cov.aiptw}{the coverage rate over \code{M} replications}
  \item{cov.tmle}{the coverage rate over \code{M} replications}
  \item{ps1}{the Min., 1st Qu., Median, Mean, 3rd Qu., Max value of probability of being treated \eqn{P(A=1|W)} for all subjects across \code{M} replications }
  \item{ps0}{the Min., 1st Qu., Median, Mean, 3rd Qu., Max value of probability of not being treated \eqn{P(A=0|W)} for all subjects across \code{M} replications }
  \item{ps11}{the Min., 1st Qu., Median, Mean, 3rd Qu., Max value of probability of being treated for subgroups with \code{A=1} across \code{M} replications }
  \item{ps00}{the Min., 1st Qu., Median, Mean, 3rd Qu., Max value of probability of not being treated for subgroups with \code{A=0}  across \code{M} replications }
}

\details{\code{print.summary.bdt}, called indirectly by entering the command \kbd{summary(result)} (where \code{result} has class \code{bdt}).}


\seealso{
\code{\link{bdt}}
}
\examples{
# Example
  set.seed(1250)
  n <- 1000
  sigma <- matrix(c(2, 1, 1, 1), ncol = 2)
  X <- matrix(rnorm(n*2), ncol = nrow(sigma)) %*% chol(sigma)
  X <- X + matrix(rep(c(.5, 1),each = n), byrow = FALSE, ncol = 2)
  I1 <- rnorm(n,mean = 1, sd = 2)
  I2 <- rnorm(n,mean = 1, sd = 1.9)
  P1 <- rnorm(n,mean = 1, sd = 1.5)
  W <- data.frame(X, I1, I2, P1)
  colnames(W) <- c("W1", "W2", "I1", "I2",  "P1")
  A <- rbinom(n, 1, plogis(0.2 + W[,"W1"] + 0.3*W[,"I1"] + W[,"W1"]*W[,"I1"]
  - 0.2*(W[,"W2"] + W[,"I2"])^2 ))
  Y <-  1 + A + W[,"W1"] + 2*W[,"W2"] + 0.5*(W[,"W1"] + W[,"P1"])^2 + rnorm(n)

  result <- bdt(Y, A, W, outcome_type = "continuous", gGLM = TRUE, gform = "A~W2+W1+I1", M = 100)
  summary(result)
}

