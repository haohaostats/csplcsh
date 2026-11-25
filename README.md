# csplcsh

**Coherent smoothed predictive likelihood for competing-risks model comparison.**

[![R-CMD-check](https://github.com/haohaostats/csplcsh/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/haohaostats/csplcsh/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE.md)

`csplcsh` compares complete cause-specific Cox models through the joint prediction
of first-event time and cause, including event-free survival.

## Installation

```r
install.packages("remotes")
remotes::install_github("haohaostats/csplcsh")
library(csplcsh)
```

Requires R >= 4.3 and the recommended package `survival`.

## Method

For each cause-specific Cox model, the package:

1. Smooths Breslow increments with a nearest-neighbour uniform kernel and boundary normalization.
2. Rescales the hazard within event-associated midpoint intervals, preserving each local Breslow increment.
3. Integrates that same adjusted hazard to obtain its cumulative hazard.

For subject $i$, the held-out joint log score is

$$
\ell_i = \sum_{j=1}^{K}\left\{
\mathbf{1}(D_i=j)\log h_j(T_i\mid x_i)-H_j(T_i\mid x_i)
\right\}.
$$

Larger mean scores are preferred. The contributions sum across causes, allowing
model differences to be traced to specific outcomes. The local increment constraint
preserves hazard mass; using the same hazard and its integral makes event densities
and event-free survival compatible.

## Quick start

Generate two competing events and independent censoring:

```r
library(csplcsh)
set.seed(2026)
n <- 180L
x <- rnorm(n)
z <- rnorm(n)
t1 <- rexp(n, 0.1 * exp(0.5 * x))
t2 <- rexp(n, 0.1 * exp(-0.5 * z))
cens <- rexp(n, 0.08)
dat <- data.frame(
  time = pmin(t1, t2, cens),
  status = ifelse(t1 < pmin(t2, cens), 1L, ifelse(t2 < cens, 2L, 0L)),
  x = x, z = z
)
```

Here `status = 0` means censoring; `1` and `2` identify the event causes. Time is
measured in years in this example. Use a consistent time unit for all candidates.

### Fit and predict

```r
fit <- cspl_fit(~ x + z, dat)
predict(fit, dat[1:3, ], times = c(0, 1, 3, 5), type = "survival")

fit_by_cause <- cspl_fit(list(`1` = ~ x, `2` = ~ z), dat)
```

`predict()` returns survival as a patient-by-time matrix. Types `"hazard"`,
`"cumhaz"`, and `"density"` return patient-by-time-by-cause arrays.
Density is $h_j(t\mid x)S(t\mid x)$; it is not cumulative incidence.

### Compare models on shared folds

```r
folds <- cspl_folds(dat$status, nfold = 5, seed = 27)
comparison <- cspl_compare(
  list(both = ~ x + z, target_only = ~ x),
  data = dat,
  folds = folds
)
comparison
comparison$cause_scores
```

The example gives approximately:

| Model | Mean CSPL–CSH | SE | Difference from leader | Paired SE |
|:--|--:|--:|--:|--:|
| both | -2.4150 | 0.0518 | 0.0000 | 0.0000 |
| target_only | -2.4517 | 0.0721 | -0.0367 | 0.0431 |

The mean cause contributions are:

| Model | Cause 1 | Cause 2 |
|:--|--:|--:|
| both | -1.2943 | -1.1207 |
| target_only | -1.2922 | -1.1595 |

In this generated example, the improvement is attributable to the competing cause.
These are tutorial results, not the simulation or application results in the paper.

For repeated cross-validation, generate common partitions internally:

```r
comparison <- cspl_compare(
  list(both = ~ x + z, target_only = ~ x), dat,
  nfold = 10, repeats = 20, seed = 2026
)
```

Alternatively pass a list of precomputed fold vectors. Every candidate uses every
partition and the same patients. A failed fit stops the comparison.

## Main functions

| Function | Purpose |
|:--|:--|
| `cspl_fit()` | Fit cause-specific Cox models and coherent baselines |
| `predict()` | Evaluate hazards, cumulative hazards, survival and event densities |
| `cspl_score()` | Return individual joint scores and cause contributions |
| `cspl_compare()` | Rank candidates with repeated, paired cross-validation |
| `cspl_cv()` | Evaluate one candidate on a specified partition |
| `cspl_folds()` | Create reproducible status-stratified folds |
| `cspl_jackknife()` | Summarize scores or paired differences |
| `cspl_truncate()` | Apply the training-endpoint rule to held-out outcomes |
| `cspl_baseline()` | Construct an adjusted baseline from Breslow increments |
| `cspl_hazard()`, `cspl_cumulative()` | Evaluate the baseline and its integral |
