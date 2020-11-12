\name{TMLE_ate}
\alias{TMLE_ate}
\title{Estimated average treatment effect by targeted maximum likelihood estimation algorithm}

\description{Using \code{tmle} function to estimate the additive effect where the outcome model is fitted by a "true" GLM and the propensity is modeled by SuperLearner (SL) and/or logistic regression. User could define the regression formula for propensity in GLM or the SL functions in SL. }

\usage{TMLE_ate(ObsData, Qform, outcome_type, gGLM, gbound = 0.025, gform = NULL, SL.library)}

\arguments{
  \item{ObsData}{a given dataset containing continuous or binary outcome \code{Y}, binary treatment \code{A}, and a set of potential confounders \code{W} based on observing n i.i.d. observations on (W,A,Y). W may contain factors. Missing values are not allowed.}
  \item{Qform}{regression formula of outcome}
  \item{gGLM}{if TRUE, GLM is used to estimate g; otherwise call \code{SuperLearner}}
  \item{gform}{optional regression formula of g; default value is "A~W"}
  \item{gbound}{value between (0,1) for truncation of predicted probabilities; default value is 0.025}
  \item{SL.library}{vector of prediction algorithms for data adaptive estimation of g}
  \item{outcome_type}{"continuous" or "binary"}
}

\value{
  \item{ATE }{a numeric value of average treatment effect computed by calling \code{tmle}}
  \item{SD }{standard error which is computed by influence function}
  \item{CI }{95\% Wald-type confidence intervals}
}
\references{
1. Gruber, S. and van der Laan, M.J. (2012), {tmle: An R Package for Targeted Maximum Likelihood Estimation}. \emph{Journal of Statistical Software}, 51(13), 1-35. \url{http://www.jstatsoft.org/v51/i13/}
}
