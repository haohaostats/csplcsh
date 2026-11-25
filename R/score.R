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

cspl_risk <- function(object, newdata, cause) {
  fit <- object$fits[[as.character(cause)]]
  if (length(stats::coef(fit)) == 0) return(rep(1, nrow(newdata)))
  r <- exp(as.numeric(stats::predict(fit, newdata = newdata, type = "lp",
                                     reference = "zero")))
  if (any(!is.finite(r))) stop("Nonfinite Cox risk multiplier")
  r
}

cspl_truncate <- function(test, tau, time = "time", status = "status") {
  .cspl_data(test, time, status)
  .cspl_scalar(tau, "tau", positive = TRUE)
  outside <- test[[time]] > tau
  test[[time]] <- pmin(test[[time]], tau)
  test[[status]][outside] <- 0L
  attr(test, "n_truncated") <- sum(outside)
  test
}

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
