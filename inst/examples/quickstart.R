library(csplcsh)

# Two competing event times and independent censoring, measured in years.
set.seed(2026)
n <- 180L
x <- rnorm(n)
z <- rnorm(n)
t1 <- rexp(n, 0.1 * exp(0.5 * x))
t2 <- rexp(n, 0.1 * exp(-0.5 * z))
cens <- rexp(n, 0.08)
dat <- data.frame(time = pmin(t1, t2, cens),
                  status = ifelse(t1 < pmin(t2, cens), 1L,
                                  ifelse(t2 < cens, 2L, 0L)), x = x, z = z)

fit <- cspl_fit(~ x + z, dat)
print(fit)
print(predict(fit, dat[1:3, ], times = c(0, 1, 3, 5)))

folds <- cspl_folds(dat$status, nfold = 5, seed = 27)
comparison <- cspl_compare(
  list(both = ~ x + z, target_only = ~ x), dat, folds = folds)
print(comparison)
print(comparison$cause_scores)

# Different predictor sets for each cause: names are the event codes.
fit_by_cause <- cspl_fit(list(`1` = ~ x, `2` = ~ z), dat)
print(fit_by_cause)
