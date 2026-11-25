\name{cspl_hazard}
\alias{cspl_hazard}
\title{Evaluate an adjusted baseline hazard}
\usage{
cspl_hazard(object, t)
}
\arguments{
\item{object}{An object from \code{\link[=cspl_baseline]{cspl_baseline()}}.}

\item{t}{Finite nonnegative evaluation times in the original time unit.}
}
\value{
A numeric vector of hazard values in the order of \code{t}.
}
\description{
Evaluate an adjusted baseline hazard
}
\details{
Internal cell boundaries belong to the cell on their right; the
training endpoint belongs to the last cell. Hazard units are inverse time.
}
\seealso{
\code{\link[=cspl_cumulative]{cspl_cumulative()}}, \code{\link[=predict.cspl_fit]{predict.cspl_fit()}}
}
