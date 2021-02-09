\name{summary.ps}
\alias{summary.ps}
\alias{print.summary.ps}
\title{Summarization of the results of a call to the ps function}
\description{These functions are \link{methods} for class \code{ps} objects, \code{summary.ps} objects}

\usage{
  \method{summary}{ps}(object, \dots)
  \method{print}{summary.ps}(x, \dots)
}

\arguments{
  \item{object}{an object of class \code{ps} for summary functions.}
  \item{x}{an object of class \code{summary.ps} for print functions.}
  \item{\dots}{currently ignored.}
}

\value{
  \item{ALL_methods}{summary of the estimated propensity scores of all methods. It is a dataframe with columns refer to Min., 1st Qu., Median, Mean, 3rd Qu., Max values and rows refer to \eqn{P(A=1|W)} and \eqn{P(A=0|W)} for all subjects, \eqn{P(A=1|W)} for subgroups \code{A=1} and \eqn{P(A=0|W)} subgroups with \code{A=0}  }
  \item{GLM1}{summary of the estimated propensity scores by user defined gform1}
  \item{GLM2}{summary of the estimated propensity scores by user defined gform2}
  \item{GLM}{summary of the estimated propensity scores if user only defined one gform, either gform1 or gform2, for the parametric estimation}
  \item{SL1}{summary of the estimated propensity scores by user defined SL.library1}
  \item{SL2}{summary of the estimated propensity scores by user defined SL.library2}
  \item{SL}{summary of the estimated propensity scores if user only defined one vector of SL.library, either SL.library1 or SL.library2}
}

\details{\code{print.summary.ps}, called indirectly by entering the command \kbd{summary(result)} (where \code{result} has class \code{ps}).}

\seealso{
  \code{\link{ps}}
  \code{\link{summary.ps}}
}
\references{
1. Bahamyirou A, Blais L, Forget A, Schnitzer ME. (2019), {Understanding and diagnosing the potential for bias when using machine learning methods with doubly robust causal estimators}. \emph{Statistical methods in medical research}, 28(6), 1637-50.
}



