# The inputs a user can produce that the mathematics cannot take, and what the
# engine currently does about them.
#
# These tests pin present behaviour, including the parts that are scheduled to
# change. When the assumption checks become returned data, this file is where
# the change becomes visible, which is the point.

test_that("uniform totals leave T undefined", {
  # If every row and column total is equal then s is that total, D has spectral
  # radius exactly 1, and I - D is singular. This is not an error by the caller.
  # A matrix of a single repeated value is the commonest way to reach it.
  A <- matrix(2, 5, 5)
  expect_warning(out <- dematel(A), "uniform-totals")
  expect_null(out)
})

test_that("a permutation matrix is degenerate for the same reason", {
  A <- matrix(c(0, 1, 0,
                0, 0, 1,
                1, 0, 0), 3, 3, byrow = TRUE)
  storage.mode(A) <- "double"
  expect_equal(rowSums(A), colSums(A))
  expect_warning(expect_null(dematel(A)), "uniform-totals")
})

test_that("the degenerate case propagates as NULL, never as an error", {
  A <- matrix(2, 4, 4)
  expect_warning(expect_null(spectral_diagnostics(A)), "uniform-totals")
  expect_warning(expect_null(sensitivity_matrix(A)),   "uniform-totals")
})

test_that("equal row totals alone are fine; it takes equal columns too", {
  # Guards the reasoning behind the fixture in test-known-answer.R. Row totals
  # all 3, column totals 3 / 2 / 4, so s = 4 and coupling is 3/4.
  A <- matrix(c(0, 2, 1,
                0, 0, 3,
                3, 0, 0), 3, 3, byrow = TRUE)
  storage.mode(A) <- "double"
  expect_true(all(rowSums(A) == 3))
  expect_false(all(colSums(A) == 3))
  expect_silent(m <- dematel(A))
  expect_equal(m$s, 4)
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
  # Two components that never meet. The Perron vector is no longer unique, so
  # the numbers mean less, but the engine must return them and let the caller
  # decide. Refusing here would make the assumption-check panel impossible.
  d <- spectral_diagnostics(block_diagonal * c(1, 2))
  expect_type(d, "list")
  expect_true(is.finite(d$lambda_max))
})

# --- Inherited behaviour that step 3 will replace ---------------------------
#
# The engine is not supposed to stop() for anything a user can cause by pasting
# an odd matrix. It currently does, because step 1 ported the five functions
# unchanged and the assumption checks are a separate step. These tests record
# the present behaviour so the replacement is a visible diff rather than a
# silent one.

test_that("INHERITED: dematel() stops on inputs a user can produce", {
  expect_error(dematel(matrix(1:6, 2, 3)))                   # not square
  expect_error(dematel(matrix(c(0, -1, 1, 0), 2, 2)))        # negative entry
  expect_error(dematel(as.data.frame(matrix(1, 2, 2))))      # not a matrix
})

test_that("INHERITED: the degenerate case signals through a warning", {
  # A warning can only be caught and re-parsed, which is how message text
  # becomes an accidental API. Step 3 turns this into a returned verdict.
  expect_warning(dematel(matrix(2, 3, 3)), "Lee et al")
})
