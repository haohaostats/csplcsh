\name{cspl_compare}
\alias{cspl_compare}
\title{Compare candidates using shared, optionally repeated cross-validation}
\usage{
cspl_compare(
  candidates,
  data,
  folds = NULL,
  nfold = 10L,
  repeats = 1L,
  seed = 1L,
  causes = NULL,
  k = 5L,
  time = "time",
  status = "status",
  order = 32
)
}
\arguments{
\item{candidates}{Named list of candidate specifications. Each entry is one
one-sided formula shared by all causes, or a list of formulas per cause.}

\item{data}{Analysis data, one row per subject.}

\item{folds}{One fold vector or a list of vectors (one per repetition).
If NULL, stratified partitions are generated using \code{seed + r - 1}.}

\item{nfold}{Number of folds for generated partitions.}

\item{repeats}{Number of generated partitions. Ignored when folds are supplied.}

\item{seed}{Nonnegative integer seed for the first generated partition.}

\item{causes}{Positive numeric event codes; NULL uses all observed causes.}

\item{k}{Number of neighbours on each side, as in \code{\link[=cspl_baseline]{cspl_baseline()}}.}

\item{time, status}{Names of follow-up and status columns.}

\item{order}{Gauss-Legendre quadrature order.}
}
\value{
A \code{cspl_comparison} object with \code{summary}, \code{cause_scores},
\code{partitions}, \code{folds}, \code{causes}, and \code{reference}. Summary rows are sorted by
decreasing mean score. Differences equal candidate minus the overall leader;
paired SEs average within-partition jackknife SEs of individual differences.
}
\description{
Compare candidates using shared, optionally repeated cross-validation
}
\details{
Larger scores are preferred. Repetitions use common folds for every
candidate. Scores are patient means; reported SEs average delete-one-fold
jackknife SEs across partitions, not the SD or SE across repetitions.
They describe variability of fold-grouped held-out scores, not full
uncertainty from training and model selection. Exact first-place ties use
original candidate order. Training-endpoint truncation is applied in each
fold. All causes must have events in every training set; errors stop the
whole comparison so every candidate uses the same observations.
}
