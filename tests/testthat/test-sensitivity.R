# The closed-form derivative, checked against finite differences.
#
# This is the check that used to live in verify_spectral_claims.R beside the
# package rather than inside it. A closed form nobody differentiates numerically
# is a closed form nobody has checked.

# lambda_max computed straight from the definition, with no shared code path
# with sensitivity_matrix(). The finite-difference check is only worth
# something if the two routes are independent.
lam_max_direct <- function(A) {
  s <- max(max(rowSums(A)), max(colSums(A)))
  D <- A / s
  T <- D %*% solve(diag(nrow(A)) - D)
  max(Re(eigen(T, only.values = TRUE)$values))
}

test_that("the closed form matches finite differences", {
  set.seed(42)
  A <- matrix(sample(0:4, 64, replace = TRUE), 8, 8)
  diag(A) <- 0
  storage.mode(A) <- "double"
  skip_if_not(is_irreducible(A))

  G <- sensitivity_matrix(A)$total

  h <- 1e-6
  FD <- matrix(NA_real_, 8, 8)
  base <- lam_max_direct(A)
  for (i in 1:8) for (j in 1:8) if (i != j) {
    Ap <- A
    Ap[i, j] <- Ap[i, j] + h
    FD[i, j] <- (lam_max_direct(Ap) - base) / h
  }

  off <- !diag(8)
  expect_lt(max(abs(G - FD)[off]), 1e-4)
})

test_that("the closed form matches finite differences across sizes", {
  set.seed(7)
  for (n in c(3, 5, 6)) {
    A <- random_admissible(n)
    G <- sensitivity_matrix(A)$total

    h <- 1e-6
    base <- lam_max_direct(A)
    for (i in seq_len(n)) for (j in seq_len(n)) if (i != j) {
      Ap <- A
      Ap[i, j] <- Ap[i, j] + h
      fd <- (lam_max_direct(Ap) - base) / h
      expect_equal(G[i, j], fd, tolerance = 1e-3,
                   info = sprintf("n = %d, entry (%d, %d)", n, i, j))
    }
  }
})

test_that("the components sum to the total and carry the right signs", {
  set.seed(11)
  for (n in c(4, 6, 9)) {
    A <- random_admissible(n)
    s <- sensitivity_matrix(A)
    off <- !diag(n)

    expect_equal(s$total[off], (s$local + s$normalization)[off], tolerance = 1e-12)

    # The local effect is always an amplification.
    expect_true(all(s$local[off] > 0))
    # The normalisation effect is never one.
    expect_true(all(s$normalization[off] <= 0))
  }
})

test_that("the normalisation effect is confined to the row or column setting s", {
  # This is what makes the counterintuitive case explicable rather than a bug:
  # a negative total sensitivity can only appear where raising the entry also
  # raises the normalising constant.
  set.seed(13)
  for (n in c(4, 6, 8)) {
    A <- random_admissible(n)
    s_const <- max(max(rowSums(A)), max(colSums(A)))
    sens <- sensitivity_matrix(A)

    critical <- matrix(FALSE, n, n)
    if (max(rowSums(A)) >= max(colSums(A)))
      critical[abs(rowSums(A) - s_const) < 1e-12, ] <- TRUE
    if (max(colSums(A)) >= max(rowSums(A)))
      critical[, abs(colSums(A) - s_const) < 1e-12] <- TRUE

    off <- !diag(n)
    expect_true(all(sens$normalization[off & !critical] == 0))
    expect_true(all(sens$total[off & !critical] > 0))
  }
})

test_that("the diagonal is NA in all three components", {
  s <- sensitivity_matrix(worked_matrices$fefo_stock_control)
  for (part in c("total", "local", "normalization")) {
    expect_true(all(is.na(diag(s[[part]]))), info = part)
    expect_false(any(is.na(s[[part]][!diag(nrow(s[[part]]))])), info = part)
  }
})

test_that("sensitivity is available for both worked matrices", {
  for (nm in names(worked_matrices)) {
    s <- sensitivity_matrix(worked_matrices[[nm]])
    expect_type(s, "list")
    expect_named(s, c("total", "local", "normalization"))
    expect_equal(dim(s$total), rep(nrow(worked_matrices[[nm]]), 2))
  }
})
