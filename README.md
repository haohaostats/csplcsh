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

## Quick start: malignant melanoma

The package includes `melanoma`, the paper-ready version of the same 205 patients
from `boot::melanoma` used in the paper:
57 melanoma deaths, 14 other-cause deaths, and 134 censored observations.
The original source is the `boot` package ([data documentation](https://stat.ethz.ch/R-manual/R-devel/library/boot/html/melanoma.html)); the processed dataset is included in `csplcsh`.

In the included dataset, follow-up has been converted from days to years and event status recoded as
0 for censoring, 1 for melanoma death, and 2 for other-cause death.
Tumour thickness is entered as `log(1 + thickness)`; the reference categories
are female sex and absent ulceration.

```r
library(csplcsh)
data("melanoma", package = "csplcsh")
dat <- melanoma
```

### Fit and predict

Model M13 includes age, log-transformed tumour thickness and ulceration.
The same predictors are fitted separately for the two causes.

```r
fit <- cspl_fit(~ age + thick + ulcer, dat)
predict(fit, dat[1:3, ], times = c(0, 1, 3, 5), type = "survival")
```

`predict()` returns survival as a patient-by-time matrix. Types `"hazard"`,
`"cumhaz"`, and `"density"` return patient-by-time-by-cause arrays.

### Compare models and decompose their scores

Compare M13 with M10, which includes thickness and ulceration. This isolates
the predictive contribution of adding age. The example uses the paper's
20 repetitions of stratified 10-fold cross-validation and the same partition seeds.

```r
candidates <- list(M13 = ~ age + thick + ulcer, M10 = ~ thick + ulcer)
comparison <- cspl_compare(
  candidates, dat,
  nfold = 10, repeats = 20, seed = 20261202
)
comparison
comparison$cause_scores
```

| Model | Predictors | Mean CSPL–CSH | SE | Difference from M13 | Paired SE |
|:--|:--|--:|--:|--:|--:|
| M13 | Age, thickness, ulceration | -1.4608 | 0.0519 | 0.0000 | 0.0000 |
| M10 | Thickness, ulceration | -1.4895 | 0.0540 | -0.0288 | 0.0159 |

The mean cause-specific contributions are:

| Model | Melanoma death | Other-cause death |
|:--|--:|--:|
| M13 | -1.0945 | -0.3663 |
| M10 | -1.0916 | -0.3979 |

Adding age improves the other-cause contribution by approximately 0.0316,
with a small decrease in the melanoma contribution. The resulting gain in
the joint score is approximately 0.0288. These values reproduce the
corresponding comparison in the paper.

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
