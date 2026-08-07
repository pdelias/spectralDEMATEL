# Analytic identities, asserted rather than assumed.
#
# Each one is a line of algebra that must hold for every admissible matrix, and
# each has been quietly wrong somewhere at some point. They run over random
# admissible matrices rather than fixtures, so they cover shapes nobody thought
# to write down.

set.seed(20260807)
admissible <- lapply(c(3, 4, 5, 7, 9, 12, 15, 20), random_admissible)

test_that("coupling, the dominant eigenvalue and the multiplier are one quantity", {
  for (A in admissible) {
    d <- spectral_diagnostics(A)
    expect_equal(d$multiplier, 1 + d$lambda_max,               tolerance = 1e-10)
    expect_equal(d$multiplier, 1 / (1 - d$mu_max),             tolerance = 1e-10)
    expect_equal(d$lambda_max, d$mu_max / (1 - d$mu_max),      tolerance = 1e-10)
    expect_gt(d$multiplier, 1)
  }
})

test_that("indirect dominance is coupling above one half", {
  for (A in admissible) {
    d <- spectral_diagnostics(A)
    expect_identical(d$indirect_dominant, d$mu_max > 0.5)
    expect_identical(d$indirect_dominant, d$lambda_max > 1)
  }
})

test_that("D and T share their eigenvectors", {
  # T = D(I - D)^-1 is a function of D, so they commute and the dominant
  # eigenvector is the same vector. The entry profile could equally have been
  # read off D.
  for (A in admissible) {
    m <- dematel(A)
    eD <- eigen(m$D)
    uD <- abs(Re(eD$vectors[, which.max(Re(eD$values))]))
    uD <- uD / sqrt(sum(uD^2))
    expect_equal(spectral_diagnostics(A)$entry_points, uD, tolerance = 1e-8)
  }
})

test_that("the two input types are independent routes to the same numbers", {
  # spectral_diagnostics() has two branches. Given A it builds T and reads the
  # spectrum; given T it recovers the spectrum of D through the inverse Moebius
  # map. Feeding it a matrix both ways exercises both branches against each
  # other, which is the cheapest independent check the package has.
  for (A in admissible) {
    from_A <- spectral_diagnostics(A, type = "A")
    from_T <- spectral_diagnostics(dematel(A)$T, type = "T")

    for (q in c("mu_max", "multiplier", "lambda_max", "dominance",
                "hierarchy_sd", "hierarchy_pr", "hierarchy_gini",
                "ev_condition")) {
      expect_equal(from_T[[q]], from_A[[q]], tolerance = 1e-8,
                   info = paste("quantity:", q))
    }
    expect_equal(from_T$entry_points, from_A$entry_points, tolerance = 1e-8)
    expect_equal(from_T$accumulation, from_A$accumulation, tolerance = 1e-8)
    expect_identical(from_T$indirect_dominant, from_A$indirect_dominant)
  }
})

test_that("every quantity stays inside its stated range", {
  for (A in admissible) {
    d <- spectral_diagnostics(A)
    n <- d$n

    expect_gt(d$mu_max, 0);          expect_lt(d$mu_max, 1)
    expect_gte(d$dominance, 0);      expect_lt(d$dominance, 1)
    expect_gt(d$hierarchy_pr, 0);    expect_lte(d$hierarchy_pr, 1 + 1e-12)
    expect_gte(d$hierarchy_gini, 0); expect_lt(d$hierarchy_gini, 1)

    # The condition number is 1 / |v'u| for two unit vectors, so it cannot be
    # below 1. The app this package replaces reported a negative number under
    # this name; that could not have survived this line.
    expect_gte(d$ev_condition, 1 - 1e-12)

    expect_equal(sqrt(sum(d$entry_points^2)), 1, tolerance = 1e-12)
    expect_equal(sqrt(sum(d$accumulation^2)), 1, tolerance = 1e-12)
    expect_length(d$entry_points, n)
    expect_length(d$accumulation, n)

    # sd of a unit vector cannot exceed the value it takes when the vector is a
    # single spike, which is the ceiling behind "hierarchy is not size-free".
    spike <- c(1, rep(0, n - 1))
    expect_lte(d$hierarchy_sd, stats::sd(spike) + 1e-12)
  }
})

test_that("mode dominance uses the modulus, not the real part", {
  # Constructed so the two disagree: the most negative eigenvalue of T carries
  # a larger modulus than the second largest real part. Reading the real part
  # gives a smaller, entirely plausible number.
  for (A in admissible) {
    lam <- eigen(dematel(A)$T)$values
    k <- which.max(Re(lam))
    by_modulus  <- max(Mod(lam[-k])) / Re(lam[k])
    by_realpart <- max(Re(lam[-k])) / Re(lam[k])

    expect_equal(spectral_diagnostics(A)$dominance, by_modulus, tolerance = 1e-10)
    if (abs(by_modulus - by_realpart) > 1e-6) {
      expect_gt(by_modulus, by_realpart)
    }
  }
})

test_that("the Gini shortcut equals the mean-absolute-difference definition", {
  # gini() uses the O(n log n) rank-weighted form. The definition is the mean
  # absolute difference over twice the mean. They must agree exactly.
  gini_by_definition <- function(x) {
    n <- length(x)
    sum(abs(outer(x, x, "-"))) / (2 * n^2 * mean(x))
  }
  for (A in admissible) {
    d <- spectral_diagnostics(A)
    expect_equal(d$hierarchy_gini, gini_by_definition(d$entry_points),
                 tolerance = 1e-10)
  }
})
