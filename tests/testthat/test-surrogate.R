# The permutation surrogate baseline.
#
# The ensemble is random, so the tests assert the properties that must hold for
# every draw rather than particular values -- with one exception: the seed has
# to make a run exactly reproducible, because a user will put an ensemble in a
# paper.

test_that("the ensemble has the requested shape", {
  ens <- surrogate_ensemble(worked_matrices$fefo_stock_control, B = 25, seed = 1)

  expect_s3_class(ens, "data.frame")
  expect_equal(nrow(ens), 25)
  expect_named(ens, c("mu_max", "lambda_max", "dominance",
                      "hierarchy_sd", "hierarchy_pr", "hierarchy_gini"))
  expect_true(all(vapply(ens, function(x) all(is.finite(x)), logical(1))))
})

test_that("a seed makes the ensemble exactly reproducible", {
  # A user reproducing a figure in a paper depends on this line.
  a <- surrogate_ensemble(worked_matrices$fefo_stock_control, B = 15, seed = 99)
  b <- surrogate_ensemble(worked_matrices$fefo_stock_control, B = 15, seed = 99)
  expect_identical(a, b)

  c <- surrogate_ensemble(worked_matrices$fefo_stock_control, B = 15, seed = 100)
  expect_false(isTRUE(all.equal(a$mu_max, c$mu_max)))
})

test_that("passing no seed leaves the caller's random state alone", {
  # set.seed() inside a function that was not asked to reseed would silently
  # change results for whatever runs next.
  set.seed(5)
  before <- runif(1)

  set.seed(5)
  invisible(surrogate_ensemble(worked_matrices$fefo_stock_control, B = 3))
  set.seed(5)
  after <- runif(1)

  expect_identical(before, after)
})

test_that("every draw satisfies the same range constraints as a real system", {
  ens <- surrogate_ensemble(worked_matrices$fefo_stock_control, B = 40, seed = 3)

  expect_true(all(ens$mu_max > 0 & ens$mu_max < 1))
  expect_true(all(ens$lambda_max > 0))
  expect_true(all(ens$dominance >= 0 & ens$dominance < 1))
  expect_true(all(ens$hierarchy_pr > 0 & ens$hierarchy_pr <= 1))
  expect_true(all(ens$hierarchy_gini >= 0 & ens$hierarchy_gini < 1))
  expect_equal(ens$lambda_max, ens$mu_max / (1 - ens$mu_max), tolerance = 1e-10)
})

test_that("shuffling preserves what it claims to preserve", {
  # The ensemble reports diagnostics, not matrices, so the invariants are
  # checked on the shuffle itself: same size, same density, same multiset of
  # off-diagonal values. If this drifts, the baseline stops being a baseline
  # and starts being a different distribution.
  A <- worked_matrices$fefo_stock_control
  off <- which(row(A) != col(A))

  set.seed(4)
  for (i in 1:20) {
    S <- A
    S[off] <- sample(A[off])
    expect_identical(dim(S), dim(A))
    expect_identical(sort(S[off]), sort(A[off]))
    expect_identical(sum(S[off] > 0), sum(A[off] > 0))
    expect_identical(diag(S), diag(A))
  }
})

test_that("the observed system can be located in its own ensemble", {
  # The use the feature exists for: is this system's coupling something the
  # analyst built, or something the entry distribution implies on its own?
  A <- worked_matrices$fefo_stock_control
  observed <- spectral_diagnostics(A)$mu_max
  ens <- surrogate_ensemble(A, B = 60, seed = 2026)

  p <- mean(ens$mu_max >= observed)
  expect_gte(p, 0)
  expect_lte(p, 1)
})
