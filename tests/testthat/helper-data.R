example_data <- function(n = 100L) {
  set.seed(101)
  x <- rnorm(n); z <- rnorm(n)
  a <- rexp(n, 0.12 * exp(0.3 * x))
  b <- rexp(n, 0.1 * exp(-0.3 * z))
  c <- rexp(n, 0.06)
  data.frame(time = pmin(a, b, c), status = ifelse(a < pmin(b, c), 1, ifelse(b < c, 2, 0)), x = x, z = z)
}
