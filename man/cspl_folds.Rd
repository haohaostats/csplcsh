\name{cspl_folds}
\alias{cspl_folds}
\title{Generate status-stratified cross-validation folds}
\usage{
cspl_folds(status, nfold = 10L, seed = 1L)
}
\arguments{
\item{status}{Numeric status vector, with 0 for censoring and positive causes.}

\item{nfold}{Integer number of nonempty folds, at least 2.}

\item{seed}{Integer seed for reproducible assignments.}
}
\value{
An integer fold vector, in the original patient order.
}
\description{
Generate status-stratified cross-validation folds
}
\details{
The random number state is restored after partition generation.
Within each status, subjects are shuffled then assigned cyclically.
Verify enough events remain in each training set; \code{\link[=cspl_cv]{cspl_cv()}} fails if a
training set has no event for any modelled cause.
}
