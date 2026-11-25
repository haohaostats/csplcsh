library(csplcsh)

data("melanoma", package = "csplcsh")
dat <- melanoma

fit <- cspl_fit(~ age + thick + ulcer, dat)
print(predict(fit, dat[1:3, ], times = c(0, 1, 3, 5), type = "survival"))

candidates <- list(M13 = ~ age + thick + ulcer, M10 = ~ thick + ulcer)
comparison <- cspl_compare(
  candidates, dat,
  nfold = 10, repeats = 20, seed = 20261202
)
print(comparison)
print(comparison$cause_scores)
