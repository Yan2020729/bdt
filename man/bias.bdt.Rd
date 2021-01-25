\name{bias}
\alias{bias.bdt}
\alias{bias}
\title{ATE bias of AIPTW and TMLE }
\description{\code{bias.bdt} is for class \code{bdt} objects which returns a dataframe of ATE bias of AIPTW and TMLE of by a call to the bdt function; \code{bias} is \link{methods} for class \code{bdt}}

\usage{
  \method{bias}{bdt}(object)
  bias(object)
}

\arguments{
  \item{object}{an object of class \code{bdt} for \code{bias} functions.}
}

\value{\item{bias.ate}{a dataframe of ATE bias with three columns: bias of AIPTW and bias of TMLE}}

\references{
1. Bahamyirou A, Blais L, Forget A, Schnitzer ME. (2019), {Understanding and diagnosing the potential for bias when using machine learning methods with doubly robust causal estimators}. \emph{Statistical methods in medical research} 28(6):1637-50.
}
\seealso{
\code{\link{bdt}}
}
