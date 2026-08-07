# The inputs a user can produce that the mathematics cannot take.
#
# The rule these pin: nothing in this package raises a condition for anything a
# user can cause by pasting an odd matrix. No stop(), no warning(), no message.
# The reason lives in assumption_checks() instead, as data.

test_that("uniform totals leave T undefined, silently", {
  # If every row and column total is equal then s is that total, D has spectral
  # radius exactly 1, and I - D is singular. Not an error by the caller: a
  # matrix in which every pair was rated identically always meets it.
  A <- matrix(2, 5, 5)
  expect_silent(out <- dematel(A))
  expect_null(out)
})

test_that("a permutation matrix is degenerate for the same reason", {
  A <- matrix(c(0, 1, 0,
                0, 0, 1,
                1, 0, 0), 3, 3, byrow = TRUE)
  storage.mode(A) <- "double"
  expect_equal(rowSums(A), colSums(A))
  expect_silent(expect_null(dematel(A)))
})

test_that("nothing a user can paste raises a condition", {
  # Every one of these used to stop() or warning(). None may now.
  inadmissible <- list(
    not_square   = matrix(1:6, 2, 3),
    negative     = matrix(c(0, -1, 1, 0), 2, 2),
    not_a_matrix = as.data.frame(matrix(1, 2, 2)),
    has_na       = matrix(c(0, NA, 1, 0), 2, 2),
    infinite     = matrix(c(0, Inf, 1, 0), 2, 2),
    all_zero     = matrix(0, 4, 4),
    uniform      = matrix(2, 4, 4),
    character    = matrix(letters[1:4], 2, 2),
    empty        = matrix(numeric(0), 0, 0)
  )

  for (nm in names(inadmissible)) {
    A <- inadmissible[[nm]]
    expect_silent(dematel(A))
    expect_null(dematel(A), info = nm)

    expect_silent(spectral_diagnostics(A))
    expect_silent(sensitivity_matrix(A))
    expect_null(sensitivity_matrix(A), info = nm)
    expect_silent(assumption_checks(A))
  }
})

test_that("an inadmissible matrix still returns the full diagnosis shape", {
  # The decision behind this: a batch job over a hundred matrices gets a
  # rectangular result and never a failure, and the numbers never travel
  # without the checks that say whether they mean anything.
  good <- spectral_diagnostics(worked_matrices$fefo_stock_control)
  bad  <- spectral_diagnostics(matrix(2, 4, 4))

  expect_identical(names(bad), names(good))
  expect_true(is.na(bad$mu_max))
  expect_true(is.na(bad$lambda_max))
  expect_true(is.na(bad$hierarchy_sd))
  expect_true(all(is.na(bad$entry_points)))
  expect_length(bad$entry_points, 4)          # n is still known
  expect_equal(bad$n, 4)

  # And it says why, in the same table shape as a passing matrix.
  expect_equal(bad$checks$verdict[bad$checks$check == "totals_vary"], "fail")
})

test_that("checks = FALSE returns NULL rather than the NA shape", {
  # The inner-loop path. surrogate_ensemble() depends on this being NULL so it
  # can reject a draw with `next`.
  expect_null(spectral_diagnostics(matrix(2, 4, 4), checks = FALSE))
  expect_null(spectral_diagnostics(matrix(1:6, 2, 3), checks = FALSE))
})

test_that("strong connectivity separates reachable from unreachable graphs", {
  expect_true(is_irreducible(cycle3))
  expect_true(is_irreducible(defective))
  expect_false(is_irreducible(block_diagonal))

  # A factor that receives but never dispatches breaks it.
  sink_node <- cycle3
  sink_node[3, ] <- 0
  expect_false(is_irreducible(sink_node))

  # As does one that dispatches but never receives.
  source_node <- cycle3
  source_node[, 1] <- 0
  expect_false(is_irreducible(source_node))

  # Only the sign pattern matters, not the weights.
  expect_identical(is_irreducible(cycle3), is_irreducible(cycle3 * 1000))
})

test_that("a disconnected graph still yields diagnostics", {
  # Two components that never meet. The entry profile is no longer unique, so
  # the numbers mean less, but the engine returns them and lets the caller
  # decide. Refusing here would make the assumption-check panel impossible.
  d <- spectral_diagnostics(block_diagonal * c(1, 2))
  expect_true(is.finite(d$lambda_max))
  expect_equal(d$checks$verdict[d$checks$check == "strong_connectivity"], "fail")
})

test_that("every diagnosis records the engine version", {
  d <- spectral_diagnostics(worked_matrices$fefo_stock_control)
  expect_identical(d$engine_version,
                   as.character(utils::packageVersion("spectralDEMATEL")))

  # Including one that could not be computed.
  expect_identical(spectral_diagnostics(matrix(2, 4, 4))$engine_version,
                   d$engine_version)
})
