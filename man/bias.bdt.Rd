\name{bias.bdt}
\alias{bias.bdt}
\alias{bias}
\title{ATE bias of IPTW, AIPTW and TMLE }
\description{\code{bias.bdt} is for class \code{bdt} objects which returns a dataframe of ATE bias of IPTW, AIPTW and TMLE of by a call to the bdt function; \code{bias} is \link{methods} for class \code{bdt}}

\usage{
  \method{bias}{bdt}(object, \dots)
  \item{\dots}{currently ignored.}
  %\method{bias}(x)
  bias(x)
}

\arguments{
  \item{object}{an object of class \code{bdt} for \code{bias.bdt} functions.}
  \item{x}{an object of class \code{bdt} for \code{bias} functions.}
  \item{\dots}{currently ignored.}
}

\value{\item{bias.ate}{a dataframe of ATE bias with three columns: bias of IPTW, bias of AIPTW and bias of TMLE}}

\references{
1. Bahamyirou A, Blais L, Forget A, Schnitzer ME. (2019), {Understanding and diagnosing the potential for bias when using machine learning methods with doubly robust causal estimators}. \emph{Statistical methods in medical research} 28(6):1637-50.
}
\seealso{
\code{\link{bdt}}
}
