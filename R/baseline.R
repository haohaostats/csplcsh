.cspl_kernel <- function(z, kernel = "uniform") {
  v <- switch(kernel,
              uniform      = 0 * z + 0.5,
              triangular   = 1 - abs(z),
              epanechnikov = 0.75 * (1 - z * z),
              quartic      = 15 / 16 * (1 - z * z)^2)
  if (is.null(v)) stop("Unknown kernel")
  v[abs(z) > 1 + 1e-14] <- 0
  pmax(v, 0)
}

.cspl_primitive <- function(z, kernel = "uniform") {
  z <- pmin(pmax(z, -1), 1)
  switch(kernel,
         uniform      = z / 2,
         triangular   = z - z * abs(z) / 2,
         epanechnikov = 0.75 * (z - z^3 / 3),
         quartic      = 15 / 16 * (z - 2 * z^3 / 3 + z^5 / 5))
}

# Raw (unscaled) k-NN kernel baseline -----------------------------------------
cspl_raw_baseline <- function(times, jumps, tau, k = 5L,
                              kernel = "uniform", tail = "error") {
  .cspl_scalar(tau, "tau", positive = TRUE)
  .cspl_scalar(k, "k", integer = TRUE, positive = TRUE)
  if (!is.numeric(times) || !is.numeric(jumps)) stop("times and jumps must be numeric")
  if (length(kernel) != 1L || is.na(kernel)) stop("Invalid kernel")
  if (length(tail) != 1L || is.na(tail)) stop("Invalid tail policy")
  if (!kernel %in% c("uniform", "triangular", "epanechnikov", "quartic"))
    stop("Unknown kernel")
  if (!tail %in% c("error", "constant")) stop("Unknown tail policy")
  if (length(k) != 1 || !is.finite(k) || k < 1 || k != as.integer(k)) stop("Invalid k")
  if (length(times) != length(jumps) || any(!is.finite(c(times, jumps, tau))) ||
      any(jumps < 0)) stop("Invalid baseline")
  keep <- jumps > 0
  times <- times[keep]; jumps <- jumps[keep]
  if (!length(times) || times[1] <= 0 || any(diff(times) <= 0) || tau < max(times))
    stop("No valid positive jumps")
  structure(list(u = times, d = jumps, tau = tau, k = as.integer(k),
                 kernel = kernel, tail = tail),
            class = "cspl_raw_baseline")
}

# k-th strictly-left / strictly-right neighbours and bandwidth b(t)
.cspl_params <- function(object, t) {
  u <- object$u; k <- object$k
  nl <- findInterval(t, u, left.open = TRUE)  # number of u < t
  nr <- findInterval(t, u)                    # number of u <= t
  li <- nl - k + 1L; ri <- nr + k
  L <- ifelse(li >= 1, u[pmax(li, 1)], 0)
  R <- ifelse(ri <= length(u), u[pmin(ri, length(u))], object$tau)
  list(L = L, R = R, b = pmax(t - L, R - t))
}

cspl_raw_hazard <- function(object, t) {
  if (any(!is.finite(t)) || any(t < 0)) stop("Invalid prediction time")
  if (object$tail == "error" && any(t > object$tau))
    stop("Prediction exceeds training follow-up")
  t <- pmin(t, object$tau); ans <- numeric(length(t))
  if (!length(t)) return(ans)
  for (start in seq.int(1, length(t), by = 1024)) {
    ids <- start:min(start + 1023, length(t)); x <- t[ids]
    b <- .cspl_params(object, x)$b
    if (any(b <= 0)) stop("Degenerate bandwidth")
    Z <- .cspl_primitive((object$tau - x) / b, object$kernel) -
         .cspl_primitive(-x / b, object$kernel)
    if (any(Z <= 0)) stop("Zero boundary mass")
    z <- outer(x, object$u, "-") / b
    ans[ids] <- as.vector(.cspl_kernel(z, object$kernel) %*% object$d) / (b * Z)
  }
  ans
}

# Break points: neighbour changes, bandwidth switches, kernel-support crossings
.cspl_breaks <- function(object, end) {
  end <- min(end, object$tau)
  edges <- sort(unique(c(0, object$u[object$u < end], end)))
  out <- edges; anchors <- c(0, object$u, object$tau)
  if (length(edges) < 2) return(edges)
  for (i in seq_len(length(edges) - 1)) {
    a <- edges[i]; b <- edges[i + 1]
    p <- .cspl_params(object, (a + b) / 2); mid <- (p$L + p$R) / 2
    if (mid > a && mid < b) out <- c(out, mid)
    hi <- min(b, mid)
    if (hi > a) { cross <- (anchors + p$R) / 2; out <- c(out, cross[cross > a & cross < hi]) }
    lo <- max(a, mid)
    if (b > lo) { cross <- (anchors + p$L) / 2; out <- c(out, cross[cross > lo & cross < b]) }
  }
  sort(unique(out))
}

.cspl_gauss_cache <- new.env(parent = emptyenv())
.cspl_gauss <- function(order) {
  .cspl_scalar(order, "order", integer = TRUE, positive = TRUE)
  key <- as.character(order)
  if (exists(key, .cspl_gauss_cache, inherits = FALSE)) return(get(key, .cspl_gauss_cache))
  if (order < 2 || order != as.integer(order)) stop("Invalid quadrature order")
  ii <- seq_len(order - 1); off <- ii / sqrt(4 * ii^2 - 1)
  J <- matrix(0, order, order); J[cbind(ii, ii + 1)] <- off; J[cbind(ii + 1, ii)] <- off
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  res <- list(x = e$values[o], w = 2 * e$vectors[1, o]^2)
  assign(key, res, .cspl_gauss_cache); res
}

cspl_raw_cumulative <- function(object, t, order = 32) {
  if (!length(t)) return(numeric())
  if (any(!is.finite(t)) || any(t < 0)) stop("Invalid prediction time")
  if (object$tail == "error" && any(t > object$tau))
    stop("Prediction exceeds training follow-up")
  x <- pmin(t, object$tau)
  if (!max(x)) return(numeric(length(x)))
  edges <- sort(unique(c(.cspl_breaks(object, max(x)), x)))
  half <- diff(edges) / 2; mid <- (head(edges, -1) + tail(edges, -1)) / 2
  q <- .cspl_gauss(order)
  points <- outer(mid, q$x, function(m, z) m) + outer(half, q$x, "*")
  hv <- matrix(cspl_raw_hazard(object, as.vector(points)), nrow = length(mid))
  areas <- half * as.vector(hv %*% q$w)
  ans <- c(0, cumsum(areas))[match(x, edges)]
  if (object$tail == "constant") ans <- ans + (t - x) * cspl_raw_hazard(object, object$tau)
  ans
}

# Locally scaled (increment-preserving) baseline -- the CSPL-CSH baseline -----
#' Construct an increment-preserving baseline hazard
#'
#' @param times Strictly increasing positive event times (zero jumps are removed).
#' @param jumps Nonnegative Breslow increments corresponding to `times`.
#' @param tau Positive training endpoint, at least the last event time.
#' @param k Positive integer number of strictly-left and strictly-right neighbours.
#'   Missing neighbours are replaced by 0 and `tau`; the bandwidth is the larger
#'   of the distances to the two selected neighbours. Default is 5.
#' @param tail `"error"` rejects evaluation beyond `tau`; `"constant"` explicitly
#'   extrapolates the terminal hazard and integrates this extension.
#' @param order Gauss-Legendre quadrature order, at least 2; default is 32.
#' @return A `cspl_baseline` object containing the raw smoother, midpoint cell
#'   boundaries, increments, local scaling factors, and integration settings.
#' @details The uniform kernel is normalized at the interval boundaries, then
#'   rescaled within event-associated midpoint cells to preserve every positive
#'   Breslow increment. The first and last boundaries are 0 and `tau`.
#'   A single positive event time is supported. Cumulative hazards integrate the
#'   same locally rescaled hazard, including its jumps at cell boundaries.
#' @examples
#' b <- cspl_baseline(c(0.8, 1.7, 2.4, 3.8, 4.6),
#'                    c(0.10, 0.13, 0.14, 0.12, 0.11), tau = 5, k = 2)
#' cspl_hazard(b, c(1, 2, 3))
#' cspl_cumulative(b, b$edges)
#' @export
cspl_baseline <- function(times, jumps, tau, k = 5L, tail = "error", order = 32) {
  base <- cspl_raw_baseline(times, jumps, tau, k = k, kernel = "uniform", tail = tail)
  u <- base$u
  edges <- c(0, (head(u, -1) + tail(u, -1)) / 2, base$tau)
  if (any(diff(edges) <= 0)) stop("Degenerate calibration cells")
  G <- cspl_raw_cumulative(base, edges, order = order)
  area <- diff(G)
  if (any(!is.finite(area)) || any(area <= 0))
    stop("Zero or invalid hazard mass in a calibration cell")
  structure(list(base = base, edges = edges, G_edges = G, d = base$d,
                 scale = base$d / area, tau = base$tau, tail = tail, order = order),
            class = "cspl_baseline")
}

.cspl_cell <- function(object, t)
  pmin(findInterval(t, object$edges, rightmost.closed = TRUE), length(object$d))

#' Evaluate an adjusted baseline hazard
#'
#' @param object An object from [cspl_baseline()].
#' @param t Finite nonnegative evaluation times in the original time unit.
#' @return A numeric vector of hazard values in the order of `t`.
#' @details Internal cell boundaries belong to the cell on their right; the
#'   training endpoint belongs to the last cell. Hazard units are inverse time.
#' @seealso [cspl_cumulative()], [predict.cspl_fit()]
#' @export
cspl_hazard <- function(object, t) {
  if (any(!is.finite(t)) || any(t < 0)) stop("Invalid time")
  if (object$tail == "error" && any(t > object$tau)) stop("Prediction beyond training range")
  x <- pmin(t, object$tau); cell <- .cspl_cell(object, x)
  cspl_raw_hazard(object$base, x) * object$scale[cell]
}

#' Evaluate the integral of an adjusted baseline hazard
#'
#' @param object An object from [cspl_baseline()].
#' @param t Finite nonnegative evaluation times in the original time unit.
#' @return Numeric cumulative hazards, with zero at time zero and cumulative
#'   sums of Breslow increments at cell boundaries.
#' @seealso [cspl_hazard()]
#' @export
cspl_cumulative <- function(object, t) {
  if (any(!is.finite(t)) || any(t < 0)) stop("Invalid time")
  if (object$tail == "error" && any(t > object$tau)) stop("Prediction beyond training range")
  x <- pmin(t, object$tau); cell <- .cspl_cell(object, x)
  local_G <- cspl_raw_cumulative(object$base, x, order = object$order) - object$G_edges[cell]
  before <- c(0, cumsum(object$d))[cell]
  H <- before + object$scale[cell] * local_G
  if (object$tail == "constant") H <- H + (t - x) * cspl_hazard(object, object$tau)
  H
}
