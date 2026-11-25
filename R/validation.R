# Validation at public boundaries; the numerical construction is unchanged.
.cspl_scalar <- function(x, name, integer = FALSE, positive = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      (integer && (x != floor(x) || abs(x) > .Machine$integer.max)) ||
      (positive && x <= 0)) stop("Invalid ", name, call. = FALSE)
}
.cspl_data <- function(data, time, status) {
  if (!is.data.frame(data) || !nrow(data)) stop("data must be a nonempty data.frame")
  if (length(time) != 1L || length(status) != 1L || anyNA(c(time, status)) ||
      time == status || !all(c(time, status) %in% names(data)))
    stop("time and status must identify distinct columns")
  if (!is.numeric(data[[time]]) || any(!is.finite(data[[time]])) || any(data[[time]] <= 0))
    stop("Follow-up times must be finite and strictly positive")
  if (!is.numeric(data[[status]]) || any(!is.finite(data[[status]])) || any(data[[status]] < 0))
    stop("Status must use 0 for censoring and positive numeric cause codes")
}
.cspl_validate_folds <- function(folds, n) {
  if (length(folds) != n || anyNA(folds) || !is.atomic(folds) ||
      !(is.numeric(folds) || is.character(folds) || is.factor(folds)) ||
      (is.numeric(folds) && any(!is.finite(folds))) || length(unique(folds)) < 2L)
    stop("folds must assign every subject to one of at least two nonempty folds")
}
