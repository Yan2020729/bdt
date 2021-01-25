\name{plot.bdt}
\alias{plot.bdt}
\title{Plotting of the results of a call to the bdt function}
\description{
  This functions is \link{methods} for class \code{bdt} objects}
\usage{
  \method{plot}{bdt}(x, xlab = NULL, ylab = "Bias", outlierSize= 0.4, notch = FALSE, pointsize = 0.7, pointcol = "blue", labelsize = 0.35, labelcol = "red", linetype = "dotted", linecol = "darkred", ...)}

\arguments{
  \item{x}{an object of class \code{bdt}.}
  \item{xlab}{x-asix label. Default value is NULL.}
  \item{ylab}{x-asix label. Default is "Bias".}
  \item{outlierSize}{the point size of outlier. Default is 0.4}
  \item{notch}{notch for boxplot. Default value is FALSE. }
  \item{pointsize}{the size of point. Default is 0.7. }
  \item{pointcol}{the color of points. Default is blue.}
  \item{labelsize}{the size of the label of coverage rates. Default is 0.35. }
  \item{labelcol}{the color of the label of coverage rates. Default is red. }
  \item{linetype}{the type of reference line. Default is dotted line. }
  \item{linecol}{the color of reference line. Default is darkred. }
  \item{...}{graphics parameters to be passed to the plotting routines.}
}

\value{
  \item{boxplot_bias_ATE}{the boxplots and points of bias of the three estimated ate and coverage rates where y-axis represents the bias of estimates, x-axis represents estimates: AIPTW, TMLE and coverage rates are marked as red in corresponding boxes}
  \item{densityplot_ps}{the 2 x 2 density plots (A, B, C, D) of the log of the estimated weights where A refers to treatment \code{A=1} in subgroups with \code{A=1}, B refers to treatment \code{A=0} in subgroups with \code{A=0}, C refers to treatment \code{A=1} for all subjects and D refers to treatment \code{A=0} for all subjects }
}
\references{
1. Bahamyirou A, Blais L, Forget A, Schnitzer ME. (2019), {Understanding and diagnosing the potential for bias when using machine learning methods with doubly robust causal estimators}. \emph{Statistical methods in medical research}, 28(6), 1637-50.
}

\seealso{
\code{\link{bdt}}
\code{\link{summary.bdt}}
}
