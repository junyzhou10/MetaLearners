% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rsim.R
\name{Rsim}
\alias{Rsim}
\title{Reference-free simplex R-learner}
\usage{
Rsim(x, y.adj, Trt, x.test, pi.hat, pi.test, W)
}
\arguments{
\item{x}{Data matrix}

\item{y.adj}{\eqn{Y - m(X)}}

\item{Trt}{Treatment assignment indicator}

\item{x.test}{Test dataset}

\item{pi.hat}{Estimated propensity score for x}

\item{pi.test}{Estimated propensity score for x.test}

\item{W}{Simplex coordinates}
}
\value{
A matrix of treatment effect mainly for optimal treatment recommendation
}
\description{
Ancillary function to obtain \code{Rsim} results
}
