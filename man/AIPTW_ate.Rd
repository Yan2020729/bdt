\name{AIPTW_ate}
\alias{AIPTW_ate}
\title{Estimated average treatment effect by augmented inverse probability treatment weighting}

\description{Using AIPTW to estimate the additive effect where the outcome model is fitted by a "true" GLM and the propensity is modeled by SuperLearner (SL) and/or logistic regression. User could define the regression formula for propensity in GLM or the SL functions in SL.
}
\usage{AIPTW_ate(ObsData, Qform, outcome_type, gGLM, gbound = 0.025, gform = NULL, SL.library)}

\arguments{
  \item{ObsData}{a given dataset containing continuous or binary outcome \code{Y}, binary treatment \code{A}, and a set of potential confounders \code{W} based on observing n i.i.d. observations on (W,A,Y). Missing values are not allowed.}
  \item{Qform}{regression formula of outcome}
  \item{gGLM}{if TRUE, GLM is used to estimate g; otherwise call \code{SuperLearner}}
  \item{gform}{optional regression formula of g; default value is "A~W"}
  \item{gbound}{value between (0,1) for truncation of predicted probabilities; default value is 0.025}
  \item{SL.library}{vector of prediction algorithms for data adaptive estimation of g}
  \item{outcome_type}{"continuous" or "binary"}
}

\value{
  \item{ATE }{a numeric value of average treatment effect}
  \item{SD }{standard error which is computed by influence function}
  \item{CI }{95\% Wald-type confidence intervals}
}

\references{
1. Bahamyirou A, Schnitzer ME, Kennedy EH, et al. Doubly robust adaptive LASSOfor effect modifier discovery. Submitted to journal.

2. van der Laan MJ and Luedtke AR. Targeted learning of an optimal dynamictreatment, and statistical inference for its  mean outcome. \emph{UC Berkeley Divisionof Biostatistics Working Paper Series}. Working Paper 2014; 317.
}
