% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/autofeature_functions.R
\name{emrselect}
\alias{emrselect}
\title{Automated feature selection}
\usage{
emrselect(dat.X, dat.S, b0 = 100, sub.boot = NULL, Scov.diag = FALSE)
}
\arguments{
\item{dat.X}{data.frame or matrix of covariates}

\item{dat.S}{data.frame or matrix of surrogate markers}

\item{b0}{number of perturbation resampling draws}

\item{Scov.diag}{if TRUE restrict covariance of S mclust to
be diagonal, default FALSE}
}
\value{
list
}
\description{
Automated feature selection
}

