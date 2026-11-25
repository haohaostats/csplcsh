test_that("local increments and the presentation example are preserved", {
  b <- cspl_baseline(c(.8, 1.7, 2.4, 3.8, 4.6), c(.1, .13, .14, .12, .11), 5, k = 2)
  expect_equal(cspl_cumulative(b, b$edges), c(0, cumsum(b$d)), tolerance = 1e-11)
  expect_equal(csplcsh:::cspl_raw_hazard(b$base, 3), .15625, tolerance = 1e-12)
  expect_true(all(cspl_hazard(b, seq(0, 5, length.out = 101)) > 0))
  expect_error(cspl_hazard(b, 6), "training")
  expect_error(cspl_baseline(1, .1, c(2, 3)), "tau")
  expect_error(cspl_baseline(1, .1, 2, order = 1), "order")
  one <- cspl_baseline(1, .2, 1)
  expect_equal(cspl_cumulative(one, c(0, 1)), c(0, .2))
  expect_length(cspl_cumulative(b, numeric()), 0L)
})
test_that("hazard integration yields normalized two-cause probability", {
  a <- cspl_baseline(c(.8, 1.7, 2.4, 3.8, 4.6), c(.1, .13, .14, .12, .11), 5, k = 2)
  b <- cspl_baseline(c(.4, 1.4, 2.8, 4.2), c(.07, .09, .1, .14), 5, k = 2)
  edges <- sort(unique(c(a$edges, b$edges,
                         csplcsh:::.cspl_breaks(a$base, 5), csplcsh:::.cspl_breaks(b$base, 5))))
  f <- function(t) (cspl_hazard(a, t) + 1.3 * cspl_hazard(b, t)) *
    exp(-cspl_cumulative(a, t) - 1.3 * cspl_cumulative(b, t))
  mass <- sum(vapply(seq_len(length(edges) - 1L), function(i)
    integrate(f, edges[i], edges[i + 1L], rel.tol = 1e-9)$value, numeric(1)))
  mass <- mass + exp(-cspl_cumulative(a, 5) - 1.3 * cspl_cumulative(b, 5))
  expect_equal(mass, 1, tolerance = 1e-8)
})
test_that("constant tail integrates the terminal hazard", {
  b <- cspl_baseline(c(1, 2, 3), c(.1, .2, .1), 4, tail = "constant")
  expect_equal(cspl_cumulative(b, 6) - cspl_cumulative(b, 4), 2 * cspl_hazard(b, 4))
})
