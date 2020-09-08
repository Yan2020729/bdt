\name{summary_ps}
\alias{summary_ps}
\title{Summarization of the propensity scores}
\description{Summarization of the probabilities of being treated and not being treated for all subjects and for subgroups \code{A=1} and \code{A=0}. Estimates could be computing by using two choices of GLM and/or SuperLearner.}

\usage{summary_ps(A, W, gform1 = NULL, gform2 = NULL, SL.library1 = NULL, SL.library2 = NULL, gbound = 0, verbose = FALSE)}

\arguments{
  \item{A}{binary treatment indicator, \code{1} - treatment, \code{0} - control}
  \item{W}{vector, matrix, or dataframe containing baseline covariates}
  \item{gform1}{optional regression formula 1 of g; default value is "A~W"}
  \item{gform2}{optional regression formula 2 of g; default value is "A~W"}
  \item{SL.library1}{vector of prediction algorithms 1 for data adaptive estimation of g}
  \item{SL.library2}{vector of prediction algorithms 2 for data adaptive estimation of g}
  \item{gbound}{bounds used to truncate g with value between (0,1) for truncation of predicted probabilities; default value is 0}
  \item{verbose}{status messages of the glm or SL summary printed if set to \code{TRUE} (default=\code{FALSE})}
}

\value{
  \item{probabilities}{a dataframe with columns refer to Min., 1st Qu., Median, Mean, 3rd Qu., Max values and rows refer to \eqn{P(A=1|W)} and \eqn{P(A=0|W)} for all subjects, \eqn{P(A=1|W)} for subgroups \code{A=1} and \eqn{P(A=0|W)} subgroups with \code{A=0}  }
  \item{summaries}{summaries of glm regression or SuperLearner}
}

\references{
1. Bahamyirou A, Blais L, Forget A, Schnitzer ME. (2019), {Understanding and diagnosing the potential for bias when using machine learning methods with doubly robust causal estimators}. \emph{Statistical methods in medical research}, 28(6), 1637-50.
}



