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

# The iteration cap. -----------------------------------------------------------
#
# The rejection loop used to run until it had B admissible draws, with no bound.
# On a sparse matrix that does not terminate in any useful time, and compiled to
# WebAssembly it is a frozen browser tab reachable by pasting a matrix.
#
# A directed cycle is the extreme: exactly as many edges as strong connectivity
# needs, so essentially no rearrangement of them survives. Measured acceptance
# is 0.000 over 400 shuffles. Adding chords buys back a small acceptance rate,
# which is what exercises the short-but-non-empty return.

cycle_matrix <- function(n = 7) {
  A <- matrix(0, n, n)
  A[cbind(seq_len(n), c(seq_len(n)[-1], 1L))] <- 1
  storage.mode(A) <- "double"
  A
}

sparse_matrix <- function() {          # measured acceptance ~0.03
  A <- cycle_matrix(7)
  A[1, 4] <- 1; A[4, 1] <- 1; A[2, 6] <- 1; A[6, 2] <- 1
  A
}

test_that("the pathological matrices really are pathological", {
  # If these stop holding, the tests below stop testing the cap and start
  # passing for the wrong reason.
  acceptance <- function(A, n = 400) {
    set.seed(1)
    off <- which(row(A) != col(A))
    mean(vapply(seq_len(n), function(i) {
      S <- A; S[off] <- sample(A[off]); is_irreducible(S)
    }, logical(1)))
  }
  expect_true(is_irreducible(cycle_matrix()))
  expect_true(is_irreducible(sparse_matrix()))
  expect_lt(acceptance(cycle_matrix()), 0.01)
  expect_lt(acceptance(sparse_matrix()), 0.20)
})

test_that("a matrix that admits no shuffle returns NULL instead of spinning", {
  elapsed <- system.time(
    ens <- surrogate_ensemble(cycle_matrix(), B = 200, seed = 1, max_attempts = 500)
  )[["elapsed"]]

  expect_null(ens)
  expect_lt(elapsed, 30)
})

test_that("a matrix that admits few shuffles returns a short ensemble that says so", {
  ens <- surrogate_ensemble(sparse_matrix(), B = 200, seed = 1, max_attempts = 400)

  expect_s3_class(ens, "data.frame")
  expect_lt(nrow(ens), 200)
  expect_gt(nrow(ens), 0)
  expect_false(attr(ens, "complete"))
  expect_equal(attr(ens, "requested"), 200)
  expect_equal(attr(ens, "attempts"), 400)
})

test_that("a complete ensemble is marked complete and does not spend its budget", {
  ens <- surrogate_ensemble(worked_matrices$fefo_stock_control, B = 20, seed = 1)

  expect_true(attr(ens, "complete"))
  expect_equal(attr(ens, "requested"), 20)
  expect_equal(nrow(ens), 20)
  expect_lt(attr(ens, "attempts"), 50L * 20L)
})

test_that("the cap does not change the draws it never truncates", {
  # A generous cap must give identical output to the uncapped original, or every
  # ensemble already published against this package silently changes.
  a <- surrogate_ensemble(worked_matrices$fefo_stock_control, B = 15, seed = 7)
  b <- surrogate_ensemble(worked_matrices$fefo_stock_control, B = 15, seed = 7,
                          max_attempts = 10000)
  expect_equal(as.data.frame(a), as.data.frame(b))
})

test_that("surrogate_position reports a short ensemble instead of hiding it", {
  sp <- surrogate_position(sparse_matrix(), B = 200, seed = 1, max_attempts = 400)

  expect_false(sp$complete)
  expect_equal(sp$requested, 200)
  expect_lt(sp$B, 200)
  expect_equal(sp$attempts, 400)
  # The share is still computed, over however many draws there were -- it is the
  # caller's job to read `complete` before quoting it.
  expect_true(is.finite(sp$type_share))
})

test_that("surrogate_position returns NULL when no shuffle is admissible", {
  expect_null(surrogate_position(cycle_matrix(), B = 200, seed = 1, max_attempts = 500))
})
