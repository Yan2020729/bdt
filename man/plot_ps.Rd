\name{plot_ps}
\alias{plot_ps}
\title{Plotting of the propensity scores}

\description{Plotting of the probabilities of being treated and not being treated for all subjects and for subgroups \code{A=1} and \code{A=0}. Estimates could be computing by using two choices of GLM and/or SuperLearner.}

\usage{plot_ps(A, W, gform1 = NULL, gform2 = NULL, SL.library1 = NULL, SL.library2 = NULL, gbound = 0,...)}

\arguments{
  \item{A}{binary treatment indicator, \code{1} - treatment, \code{0} - control}
  \item{W}{vector, matrix, or dataframe containing baseline covariates}
  \item{gform1}{optional regression formula 1 of g; default value is "A~W"}
  \item{gform2}{optional regression formula 2 of g; default value is "A~W"}
  \item{SL.library1}{vector of prediction algorithms 1 for data adaptive estimation of g}
  \item{SL.library2}{vector of prediction algorithms 2 for data adaptive estimation of g}
  \item{gbound}{value between (0,1) for truncation of predicted probabilities; default value is 0}
  \item{\dots}{currently ignored}
}

\value{
  \item{ps_plot}{the 2 x 2 density plots (A, B, C, D) of the log of the estimated weights where A refers to   treatment \code{A=1} in subgroups with \code{A=1}, B refers to treatment \code{A=0} in subgroups with \code{A=0}, C refers to treatment \code{A=1} for all subjects and D refers to treatment \code{A=0} for all subjects }
}

\references{
1. Bahamyirou A, Blais L, Forget A, Schnitzer ME. (2019), {Understanding and diagnosing the potential for bias when using machine learning methods with doubly robust causal estimators}. \emph{Statistical methods in medical research} 28(6):1637-50.
}
