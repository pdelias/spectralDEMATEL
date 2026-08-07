# Regression fixtures: the two worked matrices from the source paper's
# Section 5, pinned to nine decimal places.
#
# These exist so that a change to a definition cannot pass unnoticed. If one of
# these fails, either a definition moved on purpose -- in which case update the
# literal, bump the version, and say so in NEWS.md -- or it moved by accident,
# which is the case the file is here to catch.

test_that("FEFO stock control reproduces its published diagnostics", {
  d <- spectral_diagnostics(worked_matrices$fefo_stock_control)

  expect_equal(d$n, 13L)
  expect_equal(d$mu_max,         0.827833880044534,  tolerance = 1e-12)
  expect_equal(d$lambda_max,     4.80834371047375,   tolerance = 1e-12)
  expect_equal(d$multiplier,     5.80834371047375,   tolerance = 1e-12)
  expect_true(d$indirect_dominant)
  expect_equal(d$dominance,      0.027887448063317,  tolerance = 1e-12)
  expect_equal(d$hierarchy_sd,   0.0377567559922664, tolerance = 1e-12)
  expect_equal(d$hierarchy_pr,   0.941950487875272,  tolerance = 1e-12)
  expect_equal(d$hierarchy_gini, 0.0719461433332194, tolerance = 1e-12)
  expect_equal(d$ev_condition,   1.01465937240665,   tolerance = 1e-12)
})

test_that("Resilience capabilities reproduces its published diagnostics", {
  d <- spectral_diagnostics(worked_matrices$resilience_capabilities)

  expect_equal(d$n, 15L)
  expect_equal(d$mu_max,         0.351167023518235, tolerance = 1e-12)
  expect_equal(d$lambda_max,     0.541228692509441, tolerance = 1e-12)
  expect_equal(d$multiplier,     1.54122869250944,  tolerance = 1e-12)
  expect_false(d$indirect_dominant)
  expect_equal(d$dominance,      0.199374579781197, tolerance = 1e-12)
  expect_equal(d$hierarchy_sd,   0.147756503228298, tolerance = 1e-12)
  expect_equal(d$hierarchy_pr,   0.452239296928713, tolerance = 1e-12)
  expect_equal(d$hierarchy_gini, 0.370573454369627, tolerance = 1e-12)
  expect_equal(d$ev_condition,   2.16649708493912,  tolerance = 1e-12)
})

# The literals above and the ones below are different routes to the same three
# numbers. The values here were not computed by this package: they come from the
# precomputed spectral workbook behind the source paper's corpus, one row per
# published study, produced by separate code from separate inputs. The
# tolerances are the ones the paper's own definition guard uses.
#
# This is the check that catches a definition drifting. Agreement between a
# function and a literal it produced proves only that nothing changed;
# agreement with an independent computation proves the definition itself.

test_that("both matrices agree with the corpus, computed independently", {
  fefo <- spectral_diagnostics(worked_matrices$fefo_stock_control)
  expect_equal(fefo$lambda_max,   4.8083437105,  tolerance = 1e-4)
  expect_equal(fefo$hierarchy_sd, 0.03775675599, tolerance = 1e-5)
  expect_equal(fefo$dominance,    0.02788744806, tolerance = 1e-4)

  res <- spectral_diagnostics(worked_matrices$resilience_capabilities)
  expect_equal(res$lambda_max,   0.5412286925,  tolerance = 1e-4)
  expect_equal(res$hierarchy_sd, 0.14775650323, tolerance = 1e-5)
  expect_equal(res$dominance,    0.19937457978, tolerance = 1e-4)
})

test_that("the worked matrices straddle every threshold", {
  # If this ever fails the fixtures have been replaced with a pair that no
  # longer exercises both sides of the cuts, and the golden tests above have
  # quietly narrowed.
  fefo <- spectral_diagnostics(worked_matrices$fefo_stock_control)
  res  <- spectral_diagnostics(worked_matrices$resilience_capabilities)

  expect_gt(fefo$mu_max, 0.5)          # amplified
  expect_lt(res$mu_max,  0.5)          # dampened
  expect_lt(fefo$hierarchy_sd, 0.10)   # diffuse
  expect_gt(res$hierarchy_sd,  0.10)   # hierarchical
  expect_true(is_irreducible(worked_matrices$fefo_stock_control))
  expect_false(is_irreducible(worked_matrices$resilience_capabilities))
})

test_that("a matrix failing strong connectivity still yields diagnostics", {
  # Assumption A2 fails for this published matrix: row 14 is all zeros, so
  # factor 14 dispatches nothing. The source paper reports its diagnostics
  # anyway. The engine must not refuse -- reporting the verdict is the caller's
  # job, and refusing here would make the assumption check unimplementable.
  A <- worked_matrices$resilience_capabilities
  expect_false(is_irreducible(A))
  expect_true(all(A[14, ] == 0))

  d <- spectral_diagnostics(A)
  expect_type(d, "list")
  expect_true(is.finite(d$lambda_max))
})
