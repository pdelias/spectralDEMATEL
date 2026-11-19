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

# ---- added in 0.7.0 -----------------------------------------------------------

# Everything a user can paste, for the exports added in 0.7.0 and for the
# arguments that the older ones gained.
pasted <- function() {
  list(
    not_square   = matrix(1:6, 2, 3),
    negative     = matrix(c(0, -1, 1, 0), 2, 2),
    not_a_matrix = as.data.frame(matrix(1, 2, 2)),
    has_na       = matrix(c(0, NA, 1, 0), 2, 2),
    infinite     = matrix(c(0, Inf, 1, 0), 2, 2),
    all_zero     = matrix(0, 4, 4),
    uniform      = matrix(2, 4, 4),
    character    = matrix(letters[1:4], 2, 2),
    empty        = matrix(numeric(0), 0, 0),
    one_factor   = matrix(3, 1, 1),
    a_vector     = c(1, 2, 3),
    null         = NULL,
    disconnected = block_diagonal * c(1, 2)
  )
}

test_that("the 0.7.0 exports raise nothing for anything a user can paste", {
  for (nm in names(pasted())) {
    A <- pasted()[[nm]]
    expect_silent(normalising_line(A))
    expect_silent(prominence_relation(A))
    expect_silent(leverage_summary(A))
    expect_silent(sensitivity_matrix(A))
    for (scheme in c("full", "rows")) {
      expect_silent(surrogate_ensemble(A, B = 3, seed = 1, scheme = scheme))
      for (reading in c("H", "sd")) {
        expect_silent(surrogate_position(A, B = 3, seed = 1, scheme = scheme,
                                         reading = reading))
      }
    }
    for (reading in c("H", "sd")) {
      expect_silent(measurement_stability(A, B = 3, seed = 1, reading = reading))
      d <- spectral_diagnostics(A)
      expect_silent(structural_type(d, reading = reading))
      expect_silent(type_stability(d, reading = reading))
      expect_silent(tradeoff_residual(d, reading = reading))
    }
  }
})

test_that("the 0.7.0 exports return NULL where there is nothing to compute", {
  cannot <- c("not_square", "negative", "not_a_matrix", "has_na", "infinite",
              "all_zero", "character", "empty", "a_vector", "null")
  for (nm in cannot) {
    A <- pasted()[[nm]]
    expect_null(normalising_line(A), info = nm)
    expect_null(prominence_relation(A), info = nm)
    expect_null(leverage_summary(A), info = nm)
    expect_null(surrogate_ensemble(A, B = 3, seed = 1), info = nm)
    expect_null(surrogate_ensemble(A, B = 3, seed = 1, scheme = "rows"), info = nm)
  }

  # Uniform totals: T is undefined, so there is no prominence and no leverage.
  # The line is still a fact about the matrix.
  A <- pasted()$uniform
  expect_null(prominence_relation(A))
  expect_null(leverage_summary(A))
  expect_false(is.null(normalising_line(A)))
})

test_that("a disconnected system still has prominence, a line and a leverage summary", {
  # As with the diagnostics, the engine does not refuse. The checks say what
  # the numbers are worth.
  A <- pasted()$disconnected
  expect_s3_class(prominence_relation(A), "data.frame")
  expect_type(normalising_line(A), "list")
  expect_type(leverage_summary(A), "list")
})

test_that("an ensemble size that makes no sense returns NULL", {
  A <- worked_matrices$fefo_stock_control
  for (B in list(0, -1, NA, "ten", NULL, c(2, 3))) {
    expect_silent(out <- surrogate_ensemble(A, B = B, seed = 1))
    expect_null(out)
    expect_silent(out <- surrogate_ensemble(A, B = B, seed = 1, scheme = "rows"))
    expect_null(out)
  }
})

test_that("a cut that is not a number returns NULL", {
  d <- spectral_diagnostics(worked_matrices$fefo_stock_control)
  for (reading in c("H", "sd")) {
    expect_silent(out <- structural_type(d, hierarchy_cut = NA, reading = reading))
    expect_null(out)
    expect_silent(out <- structural_type(d, hierarchy_cut = "high", reading = reading))
    expect_null(out)
    expect_silent(out <- type_stability(d, hierarchy_cuts = character(0),
                                        reading = reading))
    expect_null(out)
  }
})
