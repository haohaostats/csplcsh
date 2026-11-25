predict.cspl_fit <- function(object, newdata, times,
                             type = c("survival", "hazard", "cumhaz", "density"), ...) {
  type <- match.arg(type)
  if (!is.data.frame(newdata) || !nrow(newdata)) stop("newdata must be nonempty")
  if (!is.numeric(times) || !length(times) || any(!is.finite(times)) || any(times < 0))
    stop("times must be finite, nonnegative and nonempty")
  dims <- c(nrow(newdata), length(times), length(object$causes))
  dn <- list(rownames(newdata), as.character(times), as.character(object$causes))
  H <- h <- array(0, dims, dimnames = dn)
  for (j in seq_along(object$causes)) {
    cc <- object$causes[j]
    r <- cspl_risk(object, newdata, cc)
    b <- object$baselines[[as.character(cc)]]
    H[, , j] <- outer(r, cspl_cumulative(b, times))
    if (type %in% c("hazard", "density")) h[, , j] <- outer(r, cspl_hazard(b, times))
  }
  S <- exp(-apply(H, c(1L, 2L), sum))
  switch(type, survival = S, hazard = h, cumhaz = H,
         density = sweep(h, c(1L, 2L), S, "*"))
}

print.cspl_fit <- function(x, ...) {
  cat("CSPL-CSH fit
")
  cat("Causes:", paste(x$causes, collapse = ", "), " | k:", x$k,
      " | training endpoint:", format(x$tau), "
")
  cat("Cox ties: Breslow; uniform kernel with local increment preservation
")
  invisible(x)
}
