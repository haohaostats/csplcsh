\name{cspl_fit}
\alias{cspl_fit}
\title{Fit coherent cause-specific Cox prediction models}
\usage{
cspl_fit(
  formulas,
  data,
  causes = NULL,
  time = "time",
  status = "status",
  k = 5L,
  tail = "error",
  order = 32
)
}
\arguments{
\item{formulas}{One one-sided formula used for all causes, or one formula
per cause. Unnamed lists follow \code{causes}; named lists use cause codes.}

\item{data}{Nonempty data frame, one row per independent subject.}

\item{causes}{Positive numeric event codes. NULL uses all observed causes.}

\item{time, status}{Names of follow-up and event-status columns. Follow-up must
be finite and strictly positive. Status 0 denotes censoring; positive codes
denote mutually exclusive causes. Missing outcomes or model predictors are
not allowed.}

\item{k}{Number of neighbours on each side; see \code{\link[=cspl_baseline]{cspl_baseline()}}.}

\item{tail}{Extrapolation policy, \code{"error"} (default) or \code{"constant"}.}

\item{order}{Gauss-Legendre quadrature order; default 32.}
}
\value{
A \code{cspl_fit} object with \code{fits} (Cox fits), \code{baselines} (adjusted
baselines), \code{raw} (Breslow cumulative hazards), \code{causes}, \code{formulas},
outcome column names, and the training endpoint \code{tau}.
}
\description{
Fit coherent cause-specific Cox prediction models
}
\details{
Fits use Breslow ties, uncentred baseline hazards, and linear
predictors referenced to zero. Each cause must have an observed event in
the training set. Cox convergence warnings should be investigated.
This interface supports ordinary right-censored, time-fixed proportional
hazards models, including interactions and transformed covariates.
Offsets, strata, clusters, frailty, delayed entry and time-varying effects
are not supported. Conditional independent censoring is assumed for the
predictive likelihood interpretation.
}
\examples{
data("melanoma", package = "csplcsh")
fit <- cspl_fit(~ age + thick + ulcer, melanoma)
predict(fit, melanoma[1:3, ], times = c(0, 1, 3, 5))
}
