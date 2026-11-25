test_that("package scores reproduce the research implementation", {
  ref <- readRDS(test_path("fixtures", "research-reference.rds"))
  d <- ref$data
  fit <- cspl_fit(~ x1 + x2, d)
  expect_equal(cspl_score(fit, d), ref$score, tolerance = 1e-10)
  actual <- cspl_cv(~ x1 + x2, d, ref$folds)
  expect_equal(actual$scores, ref$cv$scores, tolerance = 1e-10)
  expect_equal(unname(actual$components), unname(ref$cv$components), tolerance = 1e-10)
  expect_equal(actual$Lhat, ref$cv$Lhat, tolerance = 1e-10)
  expect_equal(actual$se, ref$cv$se, tolerance = 1e-10)
  for (j in seq_along(fit$causes)) {
    b <- fit$baselines[[j]]
    expect_equal(cspl_hazard(b, ref$times), ref$hazard[[j]], tolerance = 1e-10)
    expect_equal(cspl_cumulative(b, ref$times), ref$cumulative[[j]], tolerance = 1e-10)
  }
})
