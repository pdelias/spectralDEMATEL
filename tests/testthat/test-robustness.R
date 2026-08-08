# The robustness section: is the type worth anything?
#
# Both functions are stochastic, so the tests assert invariants and
# reproducibility rather than particular values -- with one exception. A seeded
# run must be exactly repeatable, because a user will put one of these numbers
# in a paper.

test_that("surrogate_position reports a position for every diagnostic", {
  sp <- surrogate_position(worked_matrices$fefo_stock_control, B = 30, seed = 1)

  expect_named(sp, c("metrics", "type", "type_share", "type_table", "B"))
  expect_equal(sp$B, 30)
  expect_setequal(sp$metrics$metric, c("mu_max", "hierarchy_sd", "dominance"))
  expect_named(sp$metrics, c("metric", "observed", "median", "min", "max",
                             "share_ge", "outside"))

  # A position is a proportion, and the summary statistics must bracket.
  expect_true(all(sp$metrics$share_ge >= 0 & sp$metrics$share_ge <= 1))
  expect_true(all(sp$metrics$min <= sp$metrics$median))
  expect_true(all(sp$metrics$median <= sp$metrics$max))

  # `outside` must agree with the range it is derived from.
  expect_identical(sp$metrics$outside,
                   sp$metrics$observed < sp$metrics$min |
                     sp$metrics$observed > sp$metrics$max)

  expect_true(sp$type_share >= 0 && sp$type_share <= 1)
  expect_equal(sum(sp$type_table), 30)
})

test_that("the observed values are the real ones, not redrawn", {
  # A surrogate summary that quietly reported an ensemble draw as the
  # observation would be invisible in every other assertion here.
  A <- worked_matrices$fefo_stock_control
  d <- spectral_diagnostics(A, checks = FALSE)
  sp <- surrogate_position(A, B = 20, seed = 7)

  for (q in c("mu_max", "hierarchy_sd", "dominance")) {
    expect_equal(sp$metrics$observed[sp$metrics$metric == q], d[[q]],
                 tolerance = 1e-12, info = q)
  }
  expect_equal(sp$type, structural_type(d)$type)
})

test_that("the FEFO system is not something its rating distribution implies", {
  # The finding the vignette reports, pinned: coupling above every draw and
  # mode dominance below every draw. If a future change to the shuffle broke
  # this, the surrogate baseline would quietly stop being a baseline.
  sp <- surrogate_position(worked_matrices$fefo_stock_control, B = 200, seed = 42)

  mu  <- sp$metrics[sp$metrics$metric == "mu_max", ]
  dom <- sp$metrics[sp$metrics$metric == "dominance", ]

  expect_true(mu$outside);  expect_equal(mu$share_ge, 0)
  expect_true(dom$outside); expect_equal(dom$share_ge, 1)

  # Hierarchy, by contrast, falls inside the ensemble on the low side.
  hier <- sp$metrics[sp$metrics$metric == "hierarchy_sd", ]
  expect_false(hier$outside)
  expect_gt(hier$share_ge, 0.5)
})

test_that("a seed makes both functions exactly reproducible", {
  A <- worked_matrices$fefo_stock_control

  a <- surrogate_position(A, B = 25, seed = 99)
  b <- surrogate_position(A, B = 25, seed = 99)
  expect_equal(a$metrics, b$metrics)
  expect_identical(a$type_share, b$type_share)

  x <- measurement_stability(A, B = 25, seed = 99)
  y <- measurement_stability(A, B = 25, seed = 99)
  expect_identical(x$share_same, y$share_same)
  expect_identical(x$type_table, y$type_table)

  # And different seeds do something different, or the seed is not being used.
  z <- measurement_stability(A, B = 25, seed = 100)
  expect_false(isTRUE(all.equal(x$type_table, z$type_table)) &&
                 isTRUE(all.equal(a$metrics$share_ge,
                                  surrogate_position(A, B = 25, seed = 100)$metrics$share_ge)))
})

test_that("measurement_stability keeps ratings inside the scale", {
  # The perturbation must not invent a negative influence, nor a rating above
  # the top of the scale the analyst used. Checked on the perturbation itself,
  # since the returned object only carries types.
  set.seed(3)
  A <- worked_matrices$fefo_stock_control
  off <- which(row(A) != col(A))
  targets <- off[A[off] > 0]

  for (i in 1:50) {
    P <- A
    P[targets] <- pmin(pmax(P[targets] + stats::runif(length(targets), -0.5, 0.5),
                            0), max(A))
    expect_true(all(P >= 0))
    expect_lte(max(P), max(A))
  }
})

test_that("recorded zeros are left alone by default and disturbed on request", {
  # A rating of zero is a judgement that there is no influence, not a one with
  # noise on it. Disturbing zeros also fills the matrix in, which destroys the
  # sparsity the connectivity checks are about.
  A <- worked_matrices$resilience_capabilities   # has genuine zeros
  expect_gt(sum(A == 0), 0)

  d <- measurement_stability(A, B = 15, seed = 5, perturb = "nonzero")
  expect_equal(d$perturb, "nonzero")

  a <- measurement_stability(A, B = 15, seed = 5, perturb = "all")
  expect_equal(a$perturb, "all")

  # Both must report which choice was made, since it changes the answer.
  expect_true(d$B_admissible > 0)
  expect_true(a$B_admissible > 0)
})

test_that("a type far from every boundary survives its own measurement noise", {
  # FEFO sits 0.33 from the coupling cut and 0.06 from the hierarchy cut, so
  # half a rating point cannot move it. If this ever fails, either the
  # perturbation grew teeth or the classification became unstable.
  ms <- measurement_stability(worked_matrices$fefo_stock_control,
                              tolerance = 0.5, B = 100, seed = 11)
  expect_equal(ms$observed_type, "diffuse-amplified")
  expect_equal(ms$share_same, 1)
  expect_equal(ms$B_admissible, 100)
})

test_that("a type sitting on the boundary is reported as the coin flip it is", {
  # This matrix has hierarchy_sd = 0.100013 against a cut of 0.10 -- a gap of
  # one part in ten thousand. Half a rating point of noise should send it either
  # way at roughly even odds, and a stability figure near 1 here would mean the
  # perturbation has no teeth.
  #
  # This replaced a monotonicity test that was both vacuous and wrong. Vacuous
  # because it compared two systems far from any boundary, where every
  # tolerance returns 1.000. Wrong because monotonicity does not hold in
  # general: a system exactly on the cut stays near 0.5 at ANY tolerance, since
  # more noise pushes it across in both directions equally.
  borderline <- matrix(c(0, 3, 0, 0, 0, 2,
                         2, 0, 1, 4, 0, 1,
                         3, 3, 0, 0, 1, 1,
                         1, 2, 3, 0, 3, 0,
                         0, 2, 3, 0, 0, 1,
                         0, 4, 3, 3, 2, 0), 6, 6, byrow = TRUE)
  storage.mode(borderline) <- "double"

  d <- spectral_diagnostics(borderline, checks = FALSE)
  expect_equal(d$hierarchy_sd, 0.10, tolerance = 1e-3)   # on the cut

  ms <- measurement_stability(borderline, tolerance = 0.5, B = 200, seed = 3)
  expect_gt(ms$share_same, 0.25)
  expect_lt(ms$share_same, 0.75)
  expect_gt(length(ms$type_table), 1)   # it really does land both sides
})

test_that("more noise does move a type that sits at a moderate distance", {
  # Monotonicity holds where there is a boundary to cross but some way to go.
  # This system's hierarchy is 0.162 against a 0.10 cut.
  core <- matrix(c(0, 4, 3, 1, 1,
                   4, 0, 3, 1, 1,
                   2, 2, 0, 1, 1,
                   1, 1, 1, 0, 1,
                   1, 1, 1, 1, 0), 5, 5, byrow = TRUE)
  storage.mode(core) <- "double"

  small <- measurement_stability(core, tolerance = 0.25, B = 200, seed = 3)
  large <- measurement_stability(core, tolerance = 1.00, B = 200, seed = 3)

  expect_equal(small$share_same, 1)
  expect_lt(large$share_same, 1)
  expect_gt(large$share_same, 0.8)
})

test_that("both return NULL rather than failing on an undiagnosable matrix", {
  degenerate <- matrix(2, 4, 4)
  expect_null(surrogate_position(degenerate, B = 5))
  expect_null(measurement_stability(degenerate, B = 5))

  # surrogate_position also declines a matrix whose surrogates cannot preserve
  # connectivity, rather than spinning in the rejection loop.
  expect_null(surrogate_position(worked_matrices$resilience_capabilities, B = 5))
})

test_that("nothing here raises a condition for a matrix a user can paste", {
  for (A in list(matrix(1:6, 2, 3), matrix(2, 4, 4),
                 matrix(c(0, -1, 1, 0), 2, 2))) {
    expect_silent(surrogate_position(A, B = 3))
    expect_silent(measurement_stability(A, B = 3))
  }
})
