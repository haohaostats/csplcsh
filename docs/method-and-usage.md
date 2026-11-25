# Method and usage guide

## Baseline construction

For event times $u_1<\cdots<u_m$ with positive Breslow increments $d_i$, let
$\tau$ be the maximum observed training follow-up. At evaluation time $t$, select
the $k$th event strictly to the left and strictly to the right. Replace unavailable
neighbours by 0 and $\tau$, respectively. If these locations are $L_k(t)$ and
$R_k(t)$, the bandwidth is $b(t)=\max\{t-L_k(t),R_k(t)-t\}$.

The raw uniform-kernel hazard is

$$
\widetilde h_0(t)=\frac{1}{b(t) Z(t)}\sum_i d_i
K\left(\frac{t-u_i}{b(t)}\right),\qquad
Z(t)=\int_{-t/b(t)}^{(\tau-t)/b(t)} K(v)\,dv,
$$

where $K(v)=\tfrac12\mathbf 1(|v|\leq1)$.
Define $e_0=0$, $e_m=\tau$, and $e_i=(u_i+u_{i+1})/2$ for interior boundaries.
For the cell associated with event $i$,

$$
a_i=\frac{d_i}{\int_{e_{i-1}}^{e_i}\widetilde h_0(s)\,ds},\qquad
h_0(t)=a_i\widetilde h_0(t),\qquad H_0(t)=\int_0^t h_0(s)\,ds.
$$

Interior boundaries belong to their right-hand cell. The endpoint belongs to the
last cell. The hazard can jump at cell boundaries; its integral is continuous.
Integration is split at neighbour changes, bandwidth switches, kernel support
crossings, cell boundaries and requested evaluation points.

For Cox relative risks $r_j(x)$, $H_j=H_{0j}r_j$, $h_j=h_{0j}r_j$,
$S=\exp(-\sum_j H_j)$ and $f_j=h_jS$. Consequently,

$$
\sum_j\int_0^\tau f_j(t\mid x)\,dt+S(\tau\mid x)=1.
$$

Local increment preservation and probability normalization are distinct properties.
Neither alone establishes optimal predictive performance for estimated models.

## Custom outcomes and cause-specific formulas

```r
fit <- cspl_fit(list(`1` = ~ x, `2` = ~ z), dat,
                time = "years", status = "event")
folds <- cspl_folds(dat$event, 5, seed = 27)
result <- cspl_compare(
  list(full = ~ x + z, separate = list(`1` = ~ x, `2` = ~ z)),
  dat, folds = folds, time = "years", status = "event"
)
```

A common formula is replicated across causes, but regression coefficients and
baselines are fitted separately. Named per-cause lists are matched by cause code.
Unnamed lists follow the supplied `causes` order, or sorted observed codes.

## External held-out scores

```r
fit <- cspl_fit(~ x + z, training)
evaluation <- cspl_truncate(validation, fit$tau)
score <- cspl_score(fit, evaluation)
mean(score$score)
colMeans(score$components)
```

`cspl_score()` does not truncate silently. `cspl_cv()` and `cspl_compare()` apply
the default endpoint rule within each training/validation split. This changes
the evaluated outcome for observations past the training endpoint, and comparisons
should use the same partitions and outcome definitions.

## Paired comparisons and uncertainty

`comparison$cause_scores` contains mean contributions by original candidate order;
their row sums equal the mean joint scores. `comparison$summary` is sorted by
score. `comparison$partitions[[r]][[m]]$scores` retains every individual score for
partition r and candidate m. Use the stored fold vectors to preserve pairing.

The jackknife deletes groups of stored held-out scores; it does not refit the
models after deleting a validation fold. For each deletion, the remaining scores
are averaged with patient weights. This accommodates unequal fold sizes.
Repeated CV reports the average score and average within-partition SE. A paired
SE is computed from patient-level candidate-minus-reference scores before
averaging over partitions. It is not the quadrature sum of separate model SEs.

## Practical scope

Rare causes need enough events for all training sets. Missing predictors and
nonfinite coefficients cause an error. Investigate any Cox convergence warnings.
An intercept-only model (`~ 1`) is supported. Standard formula transformations
and interactions are supported. Use training-defined preprocessing for externally
estimated transformations to avoid validation leakage.

Absolute log scores depend on the time unit. Comparisons require the same unit,
patients, partition and follow-up rule. The package exposes the joint criterion;
it includes the processed melanoma application data. The paper's alternative
scoring methods are not bundled.
