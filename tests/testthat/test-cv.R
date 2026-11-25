test_that("shared repeated comparison is paired and additive", {
  d <- example_data()
  folds <- list(cspl_folds(d$status, 4, 3), cspl_folds(d$status, 4, 4))
  z <- cspl_compare(list(full = ~ x + z, reduced = ~ x), d, folds = folds)
  expect_equal(rowSums(z$cause_scores)[z$summary$model], setNames(z$summary$score, z$summary$model))
  for (m in seq_len(2)) {
    leader <- match(z$reference, c("full", "reduced"))
    paired <- mean(vapply(seq_len(2), function(r)
      cspl_jackknife(z$partitions[[r]][[m]]$scores - z$partitions[[r]][[leader]]$scores, folds[[r]])$se,
      numeric(1)))
    expect_equal(z$summary$paired_se[match(c("full", "reduced")[m], z$summary$model)], paired)
  }
  expect_equal(sum(z$summary$first_percent), 100)
  tie <- cspl_compare(list(a = ~ x, b = ~ x), d, folds = folds[[1]])
  expect_equal(tie$summary$first_percent, c(100, 0))
  expect_equal(tie$summary$paired_se, c(0, 0))
})
test_that("fold generation preserves RNG and unequal-fold jackknife uses patient weights", {
  set.seed(12); before <- .Random.seed
  f <- cspl_folds(rep(0:2, 20), 5, 4)
  expect_identical(.Random.seed, before)
  expect_identical(f, cspl_folds(rep(0:2, 20), 5, 4))
  scores <- c(1, 3, 4, 8, 9, 2); folds <- c(1, 2, 2, 3, 3, 3)
  deleted <- vapply(1:3, function(r) mean(scores[folds != r]), numeric(1))
  jk <- cspl_jackknife(scores, folds)
  expect_equal(jk$Lhat, mean(scores))
  expect_equal(jk$se, sqrt(2/3 * sum((deleted - mean(deleted))^2)))
})
test_that("single-cause comparisons retain candidate rows", {
  d <- example_data(); d$status[d$status == 2] <- 0
  z <- cspl_compare(list(full = ~ x + z, reduced = ~ x), d, nfold = 4)
  expect_equal(dim(z$cause_scores), c(2L, 1L))
  expect_equal(unname(rowSums(z$cause_scores)[z$summary$model]), z$summary$score)
  one <- cspl_compare(list(full = ~ x), d, nfold = 4)
  expect_equal(dim(one$cause_scores), c(1L, 1L))
})

test_that("four causes and missing training events are handled explicitly", {
  d <- example_data(120)
  d$status <- rep(0:4, length.out = nrow(d))
  fit <- cspl_fit(~ x, d)
  expect_equal(ncol(cspl_score(fit, d)$components), 4L)
  cv <- cspl_cv(~ x, d, cspl_folds(d$status, 4))
  expect_equal(rowSums(cv$components), cv$scores)
  expect_equal(ncol(cv$components), 4L)
  d$status <- ifelse(seq_len(nrow(d)) == 1L, 2L, 1L)
  expect_error(cspl_cv(~ x, d, rep(1:4, length.out = nrow(d))), "no events")
})
