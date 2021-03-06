\name{plot.ps}
\alias{plot.ps}
\title{Plotting of the results of a call to the ps function}
\description{This functions is \link{methods} for class \code{ps} objects.}

\usage{
  \method{plot}{ps}(x, xlab_name = c("log(1/g_1)", "log(1/g_0)", "log(1/g_1)[A=1]",
  "log(1/g_0)[A=0]"), cols = c("red", "green", "blue", "orange"),...)}

\arguments{
  \item{x}{an object of class \code{ps} for plot functions.}
  \item{xlab_name}{the name of x-axis for four plots.}
  \item{cols}{the color of curves. }
  \item{...}{currently ignored.}
}

\value{
  \item{ps_plots}{the 2 x 2 density plots (A, B, C, D) of the log of the estimated weights where A refers to treatment \code{A=1} in subgroups with \code{A=1}, B refers to treatment \code{A=0} in subgroups with \code{A=0}, C refers to treatment \code{A=1} for all subjects and D refers to treatment \code{A=0} for all subjects }
}

\seealso{
  \code{\link{ps}}
  \code{\link{summary.ps}}
}

\references{
1. Bahamyirou A, Blais L, Forget A, Schnitzer ME. (2019), {Understanding and diagnosing the potential for bias when using machine learning methods with doubly robust causal estimators}. \emph{Statistical methods in medical research} 28(6):1637-50.
}
