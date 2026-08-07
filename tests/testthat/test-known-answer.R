# Matrices whose spectra can be worked out on paper.
#
# The golden files catch a definition changing. These catch a definition being
# wrong in a way that still looks plausible -- a sign, a modulus, a normalising
# constant. Each expected value below is derived independently of eigen() on the
# total-relation matrix.

test_that("pure 3-cycle: the spectrum sits on a circle, in closed form", {
  # D^3 is a multiple of the identity, because the only cycle is the whole
  # graph. Product of weights 2*3*4 = 24, normalising constant s = 4, so
  # D^3 = (24 / 4^3) I and the eigenvalues of D are the three cube roots of
  # 0.375.
  s <- 4
  expect_equal(dematel(cycle3)$s, s)

  r <- 0.375^(1 / 3)
  mu <- r * exp(2i * pi * (0:2) / 3)
  lambda <- mu / (1 - mu)
  k <- which.max(Re(lambda))

  d <- spectral_diagnostics(cycle3)
  expect_equal(d$mu_max,     Re(mu[k]),                 tolerance = 1e-12)
  expect_equal(d$lambda_max, Re(lambda[k]),             tolerance = 1e-12)
  expect_equal(d$dominance,  max(Mod(lambda[-k])) / Re(lambda[k]),
               tolerance = 1e-12)

  # The Perron vector follows from walking the cycle: v2 = 2r v1,
  # v3 = (4r/3) v2. T shares its eigenvectors with D.
  v <- c(1, 2 * r, (4 * r / 3) * (2 * r))
  expect_equal(d$entry_points, v / sqrt(sum(v^2)), tolerance = 1e-10)
})

test_that("a defective T is handled and gives the rational spectrum", {
  # Row sums 1, 1, 5 and column sums 2, 4, 1, so s = 5. The characteristic
  # polynomial of D factorises as (x - 0.4)(x + 0.2)^2, and mu = -0.2 has
  # geometric multiplicity 1, so T is not diagonalisable.
  expect_equal(dematel(defective)$s, 5)
  expect_equal(qr(dematel(defective)$D + 0.2 * diag(3))$rank, 2)

  d <- spectral_diagnostics(defective)
  expect_equal(d$mu_max,     0.4,     tolerance = 1e-12)  # mu = 0.4
  expect_equal(d$lambda_max, 2 / 3,   tolerance = 1e-12)  # 0.4 / 0.6
  expect_equal(d$multiplier, 5 / 3,   tolerance = 1e-12)  # 1 / 0.6
  expect_equal(d$dominance,  0.25,    tolerance = 1e-12)  # (1/6) / (2/3)
  expect_false(d$indirect_dominant)
})

test_that("rank-one A gives mode dominance of exactly zero", {
  # A = outer(a, b) has one non-zero eigenvalue, sum(a * b) / s, and T inherits
  # the rank. Everything below the dominant eigenvalue is zero, so the spectral
  # gap is total. This is the case where taking lambda_2 by real part rather
  # than by modulus would still look right, which is why the other fixtures
  # exist too.
  s <- 12                                  # max(row sums) = max(col sums) = 12
  expect_equal(dematel(rank1)$s, s)

  mu <- sum(rank1_a * rank1_b) / s         # 9 / 12
  expect_equal(mu, 0.75)

  d <- spectral_diagnostics(rank1)
  expect_equal(d$mu_max,     0.75,        tolerance = 1e-12)
  expect_equal(d$lambda_max, 3,           tolerance = 1e-12)  # 0.75 / 0.25
  expect_equal(d$multiplier, 4,           tolerance = 1e-12)
  expect_equal(d$dominance,  0,           tolerance = 1e-10)
  expect_true(d$indirect_dominant)

  # Entry points follow a, accumulation follows b.
  expect_equal(d$entry_points, rank1_a / sqrt(sum(rank1_a^2)), tolerance = 1e-10)
  expect_equal(d$accumulation, rank1_b / sqrt(sum(rank1_b^2)), tolerance = 1e-10)
})

test_that("2x2: every quantity in closed form", {
  # D = [[0, 1], [0.5, 0]], eigenvalues +/- sqrt(0.5).
  m <- sqrt(0.5)
  expect_equal(dematel(two_by_two)$s, 2)

  d <- spectral_diagnostics(two_by_two)
  expect_equal(d$mu_max,     m,               tolerance = 1e-12)
  expect_equal(d$lambda_max, m / (1 - m),     tolerance = 1e-12)
  # The other eigenvalue of T is -m / (1 + m), so the ratio of moduli is
  # (1 - m) / (1 + m). By real part it would have been the wrong one.
  expect_equal(d$dominance,  (1 - m) / (1 + m), tolerance = 1e-12)

  u <- c(1, m) / sqrt(1 + m^2)
  expect_equal(d$entry_points, u, tolerance = 1e-12)
  expect_equal(d$hierarchy_sd, stats::sd(u), tolerance = 1e-12)
})

test_that("hierarchy readings hit their extremes in the right direction", {
  # A perfectly even entry profile is the least concentrated system there is.
  # sd and Gini go to their minimum, the participation ratio to its maximum.
  # If a reading is ever flipped, this is where it shows.
  even <- rep(1, 6) / sqrt(6)
  expect_equal(stats::sd(even), 0)
  expect_equal((sum(even^2)^2 / sum(even^4)) / 6, 1)

  # Reaching that case needs care. Equal row totals give a uniform Perron
  # vector, but equal row AND column totals give rho(D) = 1 and no T at all
  # (see test-degenerate.R). So the row totals are held equal at 3 while a
  # column total is pushed to 4, which then sets the normalising constant.
  flat <- matrix(c(0, 2, 1,
                   0, 0, 3,
                   3, 0, 0), 3, 3, byrow = TRUE)
  storage.mode(flat) <- "double"
  expect_true(all(rowSums(flat) == 3))
  expect_equal(dematel(flat)$s, 4)

  d <- spectral_diagnostics(flat)
  expect_equal(d$mu_max, 0.75, tolerance = 1e-12)          # 3 / 4
  expect_equal(d$entry_points, rep(1, 3) / sqrt(3), tolerance = 1e-10)
  expect_equal(d$hierarchy_sd,   0, tolerance = 1e-10)   # least concentrated
  expect_equal(d$hierarchy_gini, 0, tolerance = 1e-10)   # least concentrated
  expect_equal(d$hierarchy_pr,   1, tolerance = 1e-10)   # OTHER DIRECTION

  # And the concentrated end runs the other way on all three, together.
  concentrated <- spectral_diagnostics(worked_matrices$resilience_capabilities)
  diffuse      <- spectral_diagnostics(worked_matrices$fefo_stock_control)
  expect_gt(concentrated$hierarchy_sd,   diffuse$hierarchy_sd)
  expect_gt(concentrated$hierarchy_gini, diffuse$hierarchy_gini)
  expect_lt(concentrated$hierarchy_pr,   diffuse$hierarchy_pr)   # OTHER DIRECTION
})
