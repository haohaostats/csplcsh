.cspl_formulas <- function(formulas, causes) {
  if (inherits(formulas, "formula")) formulas <- rep(list(formulas), length(causes))
  if (!is.list(formulas) || length(formulas) != length(causes))
    stop("Supply one formula or one formula per cause")
  if (!is.null(names(formulas))) {
    if (anyDuplicated(names(formulas)) || !setequal(names(formulas), as.character(causes)))
      stop("Named formulas must use the cause codes as names")
    formulas <- formulas[as.character(causes)]
  }
  for (f in formulas) {
    if (!inherits(f, "formula") || length(f) != 2L)
      stop("Use one-sided covariate formulas")
    if (any(c("offset", "strata", "cluster", "tt", "frailty") %in% all.names(f)))
      stop("Offsets, strata, clusters, frailty and time-transform terms are not supported")
  }
  names(formulas) <- as.character(causes)
  formulas
}

# Fit cause-specific Cox models and CSPL baselines on a training set
#' Fit coherent cause-specific Cox prediction models
#'
#' @param formulas One one-sided formula used for all causes, or one formula
#'   per cause. Unnamed lists follow `causes`; named lists use cause codes.
#' @param data Nonempty data frame, one row per independent subject.
#' @param causes Positive numeric event codes. NULL uses all observed causes.
#' @param time,status Names of follow-up and event-status columns. Follow-up must
#'   be finite and strictly positive. Status 0 denotes censoring; positive codes
#'   denote mutually exclusive causes. Missing outcomes or model predictors are
#'   not allowed.
#' @param k Number of neighbours on each side; see [cspl_baseline()].
#' @param tail Extrapolation policy, `"error"` (default) or `"constant"`.
#' @param order Gauss-Legendre quadrature order; default 32.
#' @return A `cspl_fit` object with `fits` (Cox fits), `baselines` (adjusted
#'   baselines), `raw` (Breslow cumulative hazards), `causes`, `formulas`,
#'   outcome column names, and the training endpoint `tau`.
#' @details Fits use Breslow ties, uncentred baseline hazards, and linear
#'   predictors referenced to zero. Each cause must have an observed event in
#'   the training set. Cox convergence warnings should be investigated.
#'   This interface supports ordinary right-censored, time-fixed proportional
#'   hazards models, including interactions and transformed covariates.
#'   Offsets, strata, clusters, frailty, delayed entry and time-varying effects
#'   are not supported. Conditional independent censoring is assumed for the
#'   predictive likelihood interpretation.
#' @examples
#' set.seed(42)
#' n <- 80
#' x <- rnorm(n)
#' t1 <- rexp(n, 0.1 * exp(0.4 * x))
#' t2 <- rexp(n, 0.1)
#' cens <- rexp(n, 0.05)
#' d <- data.frame(time = pmin(t1, t2, cens),
#'                 status = ifelse(t1 < pmin(t2, cens), 1, ifelse(t2 < cens, 2, 0)),
#'                 x = x)
#' fit <- cspl_fit(~ x, d)
#' predict(fit, d[1:2, ], times = c(0, 1, 2))
#' @export
cspl_fit <- function(formulas, data, causes = NULL, time = "time", status = "status",
                     k = 5L, tail = "error", order = 32) {
  .cspl_data(data, time, status)
  if (anyNA(data[, c(time, status)]) || any(data[[time]] <= 0))
    stop("Invalid or missing time/status")
  if (is.null(causes)) causes <- sort(unique(data[[status]][data[[status]] != 0]))
  if (!is.numeric(causes) || anyNA(causes) || anyDuplicated(causes) ||
      any(!is.finite(causes)) || !length(causes) || any(causes <= 0) || any(!data[[status]] %in% c(0, causes)))
    stop("Invalid causes")
  formulas <- .cspl_formulas(formulas, causes)
  tau <- max(data[[time]])
  fits <- baselines <- raw <- list()
  for (cc in causes) {
    id <- as.character(cc)
    dd <- data; dd$.cspl_event <- as.integer(dd[[status]] == cc)
    if (!sum(dd$.cspl_event)) stop("Training set contains no events for cause ", cc)
    f <- formulas[[id]]
    f <- stats::as.formula(call("~", as.call(list(quote(survival::Surv), as.name(time),
                              as.name(".cspl_event"))), f[[2L]]), env = environment(f))
    fit <- survival::coxph(f, data = dd, ties = "breslow", x = TRUE, y = TRUE,
                           na.action = stats::na.fail)
    if (any(!is.finite(stats::coef(fit)))) stop("Nonfinite Cox coefficient")
    bh <- survival::basehaz(fit, centered = FALSE)
    if ("strata" %in% names(bh)) stop("Stratified Cox baselines are not supported")
    inc <- c(bh$hazard[1], diff(bh$hazard))
    if (any(inc < -1e-10)) stop("Nonmonotone raw cumulative hazard")
    fits[[id]] <- fit; raw[[id]] <- bh
    baselines[[id]] <- cspl_baseline(bh$time, pmax(inc, 0), tau, k = k,
                                     tail = tail, order = order)
  }
  structure(list(fits = fits, baselines = baselines, raw = raw, causes = causes,
                 formulas = formulas, time = time, status = status, tau = tau,
                 k = as.integer(k), tail = tail),
            class = "cspl_fit")
}

# Relative risk exp(x'beta) on the zero reference, matching basehaz(centered=FALSE)
cspl_risk <- function(object, newdata, cause) {
  fit <- object$fits[[as.character(cause)]]
  if (length(stats::coef(fit)) == 0) return(rep(1, nrow(newdata)))
  r <- exp(as.numeric(stats::predict(fit, newdata = newdata, type = "lp",
                                     reference = "zero")))
  if (any(!is.finite(r))) stop("Nonfinite Cox risk multiplier")
  r
}

# Apply the tail rule to a held-out set (likelihood criteria only)
#' Apply the training-endpoint evaluation rule
#'
#' @param test Data frame containing follow-up and status columns.
#' @param tau Positive finite training endpoint.
#' @param time,status Names of follow-up and status columns.
#' @return A data frame in which follow-up beyond `tau` is truncated and marked
#'   censored (status 0). Attribute `n_truncated` counts the affected subjects.
#' @details An event exactly at `tau` retains its event status.
#' @export
cspl_truncate <- function(test, tau, time = "time", status = "status") {
  .cspl_data(test, time, status)
  .cspl_scalar(tau, "tau", positive = TRUE)
  outside <- test[[time]] > tau
  test[[time]] <- pmin(test[[time]], tau)
  test[[status]][outside] <- 0L
  attr(test, "n_truncated") <- sum(outside)
  test
}

# Per-subject held-out CSPL-CSH score, split by cause:
#   component_j = I(J = j) log h_j(T | x) - H_j(T | x);  score = sum_j component_j
#' Evaluate subject-level joint predictive log scores
#'
#' @param object A fitted object from [cspl_fit()].
#' @param test Data frame containing the fitted predictors and outcome columns.
#' @return A list with `score` (one joint log score per subject) and `components`
#'   (patient-by-cause matrix). Row sums of components equal the joint score.
#' @details Cause j contributes I(D=j) log(h_j(T|x)) - H_j(T|x).
#'   A censored observation contributes only negative cumulative hazards.
#'   Larger mean scores are preferred. This function does not truncate outcomes
#'   automatically: use [cspl_truncate()] before scoring an external test set
#'   according to the training-endpoint rule, or use [cspl_cv()].
#' @export
cspl_score <- function(object, test) {
  .cspl_data(test, object$time, object$status)
  t <- test[[object$time]]; s <- test[[object$status]]
  if (any(!s %in% c(0, object$causes))) stop("Unknown test cause")
  comp <- matrix(0, nrow(test), length(object$causes),
                 dimnames = list(NULL, paste0("cause", object$causes)))
  for (j in seq_along(object$causes)) {
    cc <- object$causes[j]; b <- object$baselines[[as.character(cc)]]
    r <- cspl_risk(object, test, cc)
    comp[, j] <- -cspl_cumulative(b, t) * r
    ev <- s == cc
    if (any(ev)) {
      h <- cspl_hazard(b, t[ev]) * r[ev]
      if (any(!(h > 0))) stop("Zero event hazard for cause ", cc)
      comp[ev, j] <- comp[ev, j] + log(h)
    }
  }
  list(score = rowSums(comp), components = comp)
}

# Status-stratified fold assignment (same rule as the DGLAS prototype)
#' Generate status-stratified cross-validation folds
#'
#' @param status Numeric status vector, with 0 for censoring and positive causes.
#' @param nfold Integer number of nonempty folds, at least 2.
#' @param seed Integer seed for reproducible assignments.
#' @return An integer fold vector, in the original patient order.
#' @details The random number state is restored after partition generation.
#'   Within each status, subjects are shuffled then assigned cyclically.
#'   Verify enough events remain in each training set; [cspl_cv()] fails if a
#'   training set has no event for any modelled cause.
#' @export
cspl_folds <- function(status, nfold = 10L, seed = 1L) {
  .cspl_scalar(nfold, "nfold", integer = TRUE, positive = TRUE)
  if (nfold < 2L || nfold > length(status)) stop("nfold must be between 2 and sample size")
  if (!is.numeric(status) || anyNA(status) || any(!is.finite(status)) || any(status < 0))
    stop("status must contain finite nonnegative numeric codes")
  .cspl_scalar(seed, "seed", integer = TRUE)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed); out <- integer(length(status))
  for (s in sort(unique(status))) {
    ids <- which(status == s); ids <- ids[sample.int(length(ids))]
    out[ids] <- rep(seq_len(nfold), length.out = length(ids))
  }
  if (length(unique(out)) != nfold) stop("Too few subjects for the requested folds")
  out
}

# Patient-weighted mean and delete-one-fold jackknife SE
#' Summarize held-out scores by a delete-one-fold jackknife
#'
#' @param scores Finite subject-level scores or paired score differences.
#' @param folds A fold vector, one entry per subject, with at least two folds.
#' @return List with patient-weighted mean `Lhat`, jackknife `se`, `fold_means`,
#'   number of folds `R`, and sample size `n`.
#' @details Deletions operate on stored held-out scores without refitting models.
#'   Unequal fold sizes retain patient weights. The SE describes variability of
#'   fold-grouped scores, not full training or model-selection uncertainty.
#' @export
cspl_jackknife <- function(scores, folds) {
  if (!is.numeric(scores) || !length(scores) || any(!is.finite(scores)))
    stop("scores must be a finite numeric vector")
  .cspl_validate_folds(folds, length(scores))
  n <- length(scores); R <- length(unique(folds))
  L <- mean(scores)
  tot <- vapply(sort(unique(folds)), function(r) sum(scores[folds == r]), numeric(1))
  nr <- as.numeric(table(factor(folds, levels = sort(unique(folds)))))
  L_minus <- (sum(scores) - tot) / (n - nr)
  var_jk <- (R - 1) / R * sum((L_minus - mean(L_minus))^2)
  list(Lhat = L, se = sqrt(var_jk), fold_means = tot / nr, R = R, n = n)
}

# Cross-validated CSPL-CSH for one candidate
#' Cross-validate a candidate using the joint predictive log score
#'
#' @param formulas Candidate specification accepted by [cspl_fit()].
#' @param data Data frame with one row per subject.
#' @param folds Fold vector with one entry per subject; at least two nonempty folds.
#' @param causes Event codes, or NULL to use all observed causes.
#' @param k Number of neighbours on each side.
#' @param tail `"error"` (default) applies training-endpoint truncation within
#'   each validation fold. `"constant"` instead uses explicit constant-hazard
#'   extrapolation and does not truncate outcomes.
#' @param time,status Names of follow-up and status columns.
#' @param order Gauss-Legendre quadrature order.
#' @return The output of [cspl_jackknife()] plus subject-level `scores`,
#'   cause-specific `components`, and total count `n_truncated`.
#' @details Supply identical folds when comparing models. A failing fold stops
#'   evaluation. [cspl_compare()] supports repeated, paired comparisons.
#' @export
cspl_cv <- function(formulas, data, folds, causes = NULL, k = 5L, tail = "error",
                    time = "time", status = "status", order = 32) {
  .cspl_data(data, time, status)
  .cspl_validate_folds(folds, nrow(data))
  if (is.null(causes)) causes <- sort(unique(data[[status]][data[[status]] != 0]))
  scores <- numeric(nrow(data))
  comp <- matrix(0, nrow(data), length(causes),
                 dimnames = list(NULL, paste0("cause", causes)))
  truncated <- 0L
  for (r in sort(unique(folds))) {
    tr <- data[folds != r, , drop = FALSE]
    obj <- cspl_fit(formulas, tr, causes = causes, k = k, tail = tail,
                    time = time, status = status, order = order)
    te <- data[folds == r, , drop = FALSE]
    if (tail == "error") { te <- cspl_truncate(te, obj$tau, time = time, status = status); truncated <- truncated + attr(te, "n_truncated") }
    sc <- cspl_score(obj, te)
    scores[folds == r] <- sc$score; comp[folds == r, ] <- sc$components
  }
  c(cspl_jackknife(scores, folds), list(scores = scores, components = comp,
                                        n_truncated = truncated))
}
