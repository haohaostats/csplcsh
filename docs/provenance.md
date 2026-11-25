# Numerical provenance

The version 0.1.0 package is derived from the paper's analysis implementation.
The reference fixture uses synthetic data only; it contains no patient records.

## Original source fingerprints (MD5)

- `cspl_smoother.R`: `be1406dc5704d477d4e2f4617ffc4e58`
- `cspl_score.R`: `c445841026380d74134d70d2a04e724e`
- `dgp.R`: `2b72c6ae775ae743d4f214854721cc7b`

## Fixture generation

R 4.5.2; set.seed(39); generate_nonnested(120, 0.6). Follow-up times
were rounded to one decimal and bounded below at 0.1 to exercise ties.
cspl_folds(status, 5, 913); cspl_fit(~ x1 + x2); default k = 5 and order = 32.
Stored outputs are individual scores and cause components, five-fold CV
scores and jackknife summaries, and hazards/integrals on 37 evaluation times.

## Package-specific interface additions

Input checks, named cause-formula matching, custom outcome column names,
RNG preservation for fold generation, prediction arrays, and repeated paired
comparison are package interfaces. The smoothing, local scaling, integration,
Breslow fitting convention, likelihood and jackknife calculations follow the research code.
No simulation or application result files were modified when creating this package.
