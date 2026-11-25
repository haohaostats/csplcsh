#' Compare candidates using shared, optionally repeated cross-validation
#'
#' @param candidates Named list of candidate specifications. Each entry is one
#'   one-sided formula shared by all causes, or a list of formulas per cause.
#' @param data Analysis data, one row per subject.
#' @param folds One fold vector or a list of vectors (one per repetition).
#'   If NULL, stratified partitions are generated using `seed + r - 1`.
#' @param nfold Number of folds for generated partitions.
#' @param repeats Number of generated partitions. Ignored when folds are supplied.
#' @param seed Nonnegative integer seed for the first generated partition.
#' @param causes Positive numeric event codes; NULL uses all observed causes.
#' @param k Number of neighbours on each side, as in [cspl_baseline()].
#' @param time,status Names of follow-up and status columns.
#' @param order Gauss-Legendre quadrature order.
#' @return A `cspl_comparison` object with `summary`, `cause_scores`,
#'   `partitions`, `folds`, `causes`, and `reference`. Summary rows are sorted by
#'   decreasing mean score. Differences equal candidate minus the overall leader;
#'   paired SEs average within-partition jackknife SEs of individual differences.
#' @details Larger scores are preferred. Repetitions use common folds for every
#'   candidate. Scores are patient means; reported SEs average delete-one-fold
#'   jackknife SEs across partitions, not the SD or SE across repetitions.
#'   They describe variability of fold-grouped held-out scores, not full
#'   uncertainty from training and model selection. Exact first-place ties use
#'   original candidate order. Training-endpoint truncation is applied in each
#'   fold. All causes must have events in every training set; errors stop the
#'   whole comparison so every candidate uses the same observations.
#' @export
cspl_compare <- function(candidates, data, folds = NULL, nfold = 10L,
                         repeats = 1L, seed = 1L, causes = NULL, k = 5L,
                         time = "time", status = "status", order = 32) {
  .cspl_data(data, time, status)
  if (!is.list(candidates) || !length(candidates)) stop("Supply a list of candidates")
  nm <- names(candidates)
  if (is.null(nm)) nm <- paste0("M", seq_along(candidates))
  if (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm)) stop("Candidate names must be unique")
  names(candidates) <- nm
  if (is.null(causes)) causes <- sort(unique(data[[status]][data[[status]] != 0]))
  if (is.null(folds)) {
    .cspl_scalar(repeats, "repeats", integer = TRUE, positive = TRUE)
    folds <- lapply(seq_len(repeats), function(r) cspl_folds(data[[status]], nfold, seed + r - 1L))
  } else if (!is.list(folds)) folds <- list(folds)
  if (!length(folds)) stop("Supply at least one partition")
  for (f in folds) .cspl_validate_folds(f, nrow(data))
  partitions <- lapply(folds, function(f) lapply(candidates, function(spec)
    cspl_cv(spec, data, f, causes, k, time = time, status = status, order = order)))
  means <- vapply(partitions, function(p) vapply(p, `[[`, numeric(1), "Lhat"), numeric(length(nm)))
  means <- matrix(means, nrow = length(nm), dimnames = list(nm, NULL))
  overall <- rowMeans(means)
  best <- which.max(overall)
  se <- vapply(seq_along(nm), function(m)
    mean(vapply(partitions, function(p) p[[m]]$se, numeric(1))), numeric(1))
  paired <- vapply(seq_along(nm), function(m)
    mean(vapply(seq_along(partitions), function(r) {
      p <- partitions[[r]]
      cspl_jackknife(p[[m]]$scores - p[[best]]$scores, folds[[r]])$se
    }, numeric(1))), numeric(1))
  first <- tabulate(apply(means, 2L, which.max), nbins = length(nm)) / length(folds) * 100
  summary <- data.frame(model = nm, score = overall, se = se,
                        difference = overall - overall[best], paired_se = paired,
                        first_percent = first, row.names = NULL)
  summary <- summary[order(-summary$score, seq_along(nm)), , drop = FALSE]
  rownames(summary) <- NULL
  cause_scores <- vapply(seq_along(nm), function(m)
    Reduce(`+`, lapply(partitions, function(p) colMeans(p[[m]]$components))) / length(folds),
    numeric(length(causes)))
  cause_scores <- t(matrix(cause_scores, nrow = length(causes), ncol = length(nm)))
  dimnames(cause_scores) <- list(nm, paste0("cause", causes))
  structure(list(summary = summary, cause_scores = cause_scores,
                 partitions = partitions, folds = folds, causes = causes,
                 reference = nm[best]), class = "cspl_comparison")
}

#' @export
print.cspl_comparison <- function(x, ...) {
  cat("CSPL-CSH model comparison (higher scores preferred)
")
  print(x$summary, row.names = FALSE)
  cat("Reference:", x$reference, " | partitions:", length(x$partitions), "
")
  invisible(x)
}
