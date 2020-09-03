\name{plot.bdt}
\alias{plot.bdt}
\title{Plotting of the results of a call to the bdt function}
\description{
  This functions is \link{methods} for class \code{bdt} objects}
\usage{
  \method{plot}{bdt}(object, \dots)}

\arguments{
  \item{object}{an object of class \code{bdt}.}
  \item{\dots}{currently ignored.}
}

\value{
  \item{boxplot_bias_ATE}{the boxplots and points of bias of the three estimated ate and coverage rates where y-axis represents the bias of estimates, x-axis represents estimates: IPTW, AIPTW, TMLE and coverage rates are marked as red in corresponding boxes}
  \item{densityplot_ps}{the 2 x 2 density plots (A, B, C, D) of the log of the estimated weights where A refers to treatment \code{A=1} in subgroups with \code{A=1}, B refers to treatment \code{A=0} in subgroups with \code{A=0}, C refers to treatment \code{A=1} for all subjects and D refers to treatment \code{A=0} for all subjects }
}
\references{
1. Bahamyirou A, Blais L, Forget A, Schnitzer ME. (2019), {Understanding and diagnosing the potential for bias when using machine learning methods with doubly robust causal estimators}. \emph{Statistical methods in medical research}, 28(6), 1637-50.
}

\seealso{
\code{\link{bdt}}

\code{\link{summary.bdt}}
}
