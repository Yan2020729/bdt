\name{summary.ps}
\alias{summary.ps}
\title{Summarization of the results of a call to the ps function}
\description{These functions are \link{methods} for class \code{ps} objects}

\usage{
  \method{summary}{ps}(object, \dots)}

\arguments{
  \item{object}{an object of class \code{ps} for summary functions.}
  \item{\dots}{currently ignored.}
}

\value{
  \item{summ_list}{a dataframe with columns refer to Min., 1st Qu., Median, Mean, 3rd Qu., Max values and rows refer to \eqn{P(A=1|W)} and \eqn{P(A=0|W)} for all subjects, \eqn{P(A=1|W)} for subgroups \code{A=1} and \eqn{P(A=0|W)} subgroups with \code{A=0}  }
}
\seealso{
  \code{\link{ps}}
  \code{\link{summary.ps}}
}
\references{
1. Bahamyirou A, Blais L, Forget A, Schnitzer ME. (2019), {Understanding and diagnosing the potential for bias when using machine learning methods with doubly robust causal estimators}. \emph{Statistical methods in medical research}, 28(6), 1637-50.
}



