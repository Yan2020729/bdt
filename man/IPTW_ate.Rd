\name{IPTW_ate}
\alias{IPTW_ate}
\title{Estimated average treatment effect by inverse probability treatment weighting}

\description{Using IPTW to estimate the additive effect where the outcome model is fitted by a "true" GLM and the propensity is modeled by SuperLearner (SL) and/or logistic regression. User could define the regression formula for propensity in GLM or the SL functions in SL.
}
\usage{IPTW_ate(ObsData, outcome_type, gGLM, gbound = 0.025, gform = NULL, SL.library, gMeth = FALSE)}

\arguments{
  \item{ObsData}{a given dataset containing continuous or binary outcome \code{Y}, binary treatment \code{A}, and a set of potential confounders \code{W} based on observing n i.i.d. observations on (W,A,Y). Missing values are not allowed.}
  \item{gGLM}{if TRUE, GLM is used to estimate g; otherwise call \code{SuperLearner}}
  \item{gform}{optional regression formula of g; default value is "A~W"}
  \item{gbound}{value between (0,1) for truncation of predicted probabilities; default value is 0.025}
  \item{SL.library}{vector of prediction algorithms for data adaptive estimation of g}
  \item{outcome_type}{"continuous" or "binary"}
  \item{gMeth}{if TRUE, print status message of the estimation method for g; default value is FALSE}
}

\value{
  \item{ATE }{a numeric value of average treatment effect}
  \item{SD}{standard error which is computed by calling \code{sandwich}}
  \item{CI }{95\% confidence intervals computed}
}

\references{
1. Zeileis A (2004). Econometric computing with HC and HAC covariance matrix estimators. \emph{Journal of Statistical Software}, 11(10), 1-17.
}
