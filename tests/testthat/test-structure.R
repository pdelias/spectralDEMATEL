# The structure map.
#
# The corpus constants are the risky part: they were extracted once, from a
# pipeline that no longer runs here, and nothing else in the package would
# notice if one were mistyped. Several are therefore checked against each other
# and against figures the source paper prints, which is an independent route.

test_that("the corpus constants are internally consistent", {
  cp <- spectralDEMATEL:::CORPUS

  expect_equal(sum(cp$count), cp$n)
  expect_equal(cp$n, 117L)
  expect_setequal(names(cp$count), names(cp$multiplier))
  expect_setequal(names(cp$count), names(spectralDEMATEL:::INTERVENTION_LOGIC))

  # Shares the source paper quotes: 69.2 / 13.7 / 10.3 / 6.8 per cent.
  share <- round(100 * cp$count / cp$n, 1)
  expect_equal(unname(share["diffuse-amplified"]),      69.2)
  expect_equal(unname(share["hierarchical-dampened"]),  13.7)
  expect_equal(unname(share["hierarchical-amplified"]), 10.3)
  expect_equal(unname(share["diffuse-dampened"]),        6.8)

  # Quantiles must be sorted, and the median must sit between the quartiles.
  for (q in cp$quantiles) {
    expect_identical(q, sort(q))
    expect_gt(q[3], q[2]); expect_lt(q[3], q[4])
  }
  expect_length(cp$probs, 5)

  # The trade-off runs downward and explains what the paper says it explains.
  expect_lt(cp$tradeoff$slope, 0)
  expect_gt(cp$tradeoff$resid_sd, 0)
  expect_equal(cp$tradeoff$r2, 0.591, tolerance = 1e-3)
})

test_that("the fitted line passes through the corpus medians, near enough", {
  # An independent route to the same constants: a least-squares line through a
  # cloud passes through its centroid. The medians are not the means, so this
  # is a sanity bound rather than an identity -- but a mistyped intercept or a
  # sign error on the slope would miss by far more than this.
  cp <- spectralDEMATEL:::CORPUS
  med_mu   <- cp$quantiles$mu_max[3]
  med_hier <- cp$quantiles$hierarchy_sd[3]
  predicted <- cp$tradeoff$intercept + cp$tradeoff$slope * med_mu
  expect_lt(abs(predicted - med_hier), 0.5 * cp$tradeoff$resid_sd)
})

test_that("the four types are named from the two cuts", {
  mk <- function(mu, hier) list(mu_max = mu, hierarchy_sd = hier)

  expect_equal(structural_type(mk(0.80, 0.05))$type, "diffuse-amplified")
  expect_equal(structural_type(mk(0.80, 0.20))$type, "hierarchical-amplified")
  expect_equal(structural_type(mk(0.30, 0.20))$type, "hierarchical-dampened")
  expect_equal(structural_type(mk(0.30, 0.05))$type, "diffuse-dampened")

  # Exactly on a cut falls to the lower side, and the margin is zero.
  s <- structural_type(mk(0.5, 0.10))
  expect_equal(s$type, "diffuse-dampened")
  expect_equal(s$coupling_margin, 0)
  expect_equal(s$hierarchy_margin, 0)
})

test_that("the margin is signed and points the right way", {
  mk <- function(mu, hier) list(mu_max = mu, hierarchy_sd = hier)

  s <- structural_type(mk(0.95, 0.18))
  expect_gt(s$coupling_margin, 0)     # amplified
  expect_gt(s$hierarchy_margin, 0)    # hierarchical

  s <- structural_type(mk(0.20, 0.02))
  expect_lt(s$coupling_margin, 0)     # dampened
  expect_lt(s$hierarchy_margin, 0)    # diffuse

  # The whole point of the margin: two systems of the same type, one of them a
  # far closer call than the other.
  near <- structural_type(mk(0.51, 0.05))
  far  <- structural_type(mk(0.95, 0.05))
  expect_equal(near$type, far$type)
  expect_lt(near$nearest_margin_scaled, far$nearest_margin_scaled)
})

test_that("nearest_cut compares the axes on a common scale", {
  # Raw distances would always name coupling, because coupling spans (0,1) while
  # hierarchy_sd rarely passes 0.2. Scaling by each axis's corpus interquartile
  # range is what makes the comparison mean anything.
  mk <- function(mu, hier) list(mu_max = mu, hierarchy_sd = hier)

  # Far from the coupling cut, a whisker from the hierarchy cut.
  s <- structural_type(mk(0.95, 0.101))
  expect_equal(s$nearest_cut, "hierarchy")
  expect_lt(s$nearest_margin_scaled, 0.05)

  # And the other way round.
  s <- structural_type(mk(0.505, 0.01))
  expect_equal(s$nearest_cut, "coupling")

  # Raw distance would have got the first case wrong: 0.001 in hierarchy is a
  # smaller number than 0.45 in coupling, but only the scaled comparison knows
  # that 0.001 is also a smaller share of hierarchy's spread.
  expect_lt(abs(0.101 - 0.10), abs(0.95 - 0.5))
})

test_that("the worked matrices land where the corpus says they do", {
  fefo <- structural_type(spectral_diagnostics(worked_matrices$fefo_stock_control))
  expect_equal(fefo$type, "diffuse-amplified")
  expect_equal(fefo$corpus_share, 81 / 117)

  resl <- structural_type(spectral_diagnostics(worked_matrices$resilience_capabilities))
  expect_equal(resl$type, "hierarchical-dampened")
  expect_equal(resl$corpus_share, 16 / 117)

  # Between them they occupy opposite corners, which is why they are the
  # worked pair.
  expect_gt(fefo$coupling_margin, 0);  expect_lt(fefo$hierarchy_margin, 0)
  expect_lt(resl$coupling_margin, 0);  expect_gt(resl$hierarchy_margin, 0)
})

test_that("the intervention logic is quoted, and never travels without its caveat", {
  for (ty in names(spectralDEMATEL:::INTERVENTION_LOGIC)) {
    txt <- spectralDEMATEL:::INTERVENTION_LOGIC[[ty]]
    expect_gt(nchar(txt), 60)
    # It describes an appraisal basis and a spread of effort, never an
    # imperative addressed to the reader.
    expect_false(grepl("\\byou\\b|\\byour\\b|\\bmust\\b|\\bshould do\\b", txt))
  }

  s <- structural_type(spectral_diagnostics(worked_matrices$fefo_stock_control))
  expect_true(nzchar(s$caveat))
  expect_match(s$caveat, "hypothesis")
  expect_match(s$caveat, "not been\\s+validated")
})

test_that("both cuts are reported, and whether each was the default", {
  d <- spectral_diagnostics(worked_matrices$fefo_stock_control)

  s <- structural_type(d)
  expect_true(s$cuts$coupling_is_default)
  expect_true(s$cuts$hierarchy_is_default)

  # The user-defined cut is a settled feature, so moving it must be reported
  # rather than silently absorbed -- and must actually change the verdict when
  # it crosses the value. FEFO's hierarchy_sd is 0.0378.
  s2 <- structural_type(d, hierarchy_cut = 0.05)     # still above the value
  expect_false(s2$cuts$hierarchy_is_default)
  expect_equal(s2$cuts$hierarchy, 0.05)
  expect_equal(s2$type, "diffuse-amplified")

  s3 <- structural_type(d, hierarchy_cut = 0.03)     # now below it
  expect_equal(s3$type, "hierarchical-amplified")
})

test_that("classification always states which hierarchy reading it used", {
  # Three readings, two directions. A type carries no meaning without this.
  s <- structural_type(spectral_diagnostics(worked_matrices$fefo_stock_control))
  expect_equal(s$hierarchy_reading, "hierarchy_sd")
})

test_that("type_stability answers the question the cut being a recommendation raises", {
  # FEFO is far from the hierarchy cut (0.038 against 0.10), so no cut in the
  # recommended range can move it.
  st <- type_stability(spectral_diagnostics(worked_matrices$fefo_stock_control))
  expect_true(st$stable)
  expect_true(all(st$by_cut$type == "diffuse-amplified"))
  expect_true(all(is.na(st$flips_between)))

  # A system placed deliberately inside the range must flip, and must say
  # between which two cuts.
  borderline <- list(mu_max = 0.8, hierarchy_sd = 0.095)
  st2 <- type_stability(borderline)
  expect_false(st2$stable)
  expect_length(st2$flips_between, 2)
  expect_lt(st2$flips_between[1], 0.095)
  expect_gt(st2$flips_between[2], 0.095)
})

test_that("the trade-off residual is signed and scaled", {
  d <- spectral_diagnostics(worked_matrices$resilience_capabilities)
  r <- tradeoff_residual(d)

  cp <- spectralDEMATEL:::CORPUS
  expect_equal(r$expected, cp$tradeoff$intercept + cp$tradeoff$slope * d$mu_max)
  expect_equal(r$residual, d$hierarchy_sd - r$expected)
  expect_equal(r$residual_sd, r$residual / cp$tradeoff$resid_sd)
  expect_equal(r$direction, if (r$residual >= 0) "above" else "below")

  # And the finding worth pinning: this system is one of the two UNCOMMON
  # types, yet it sits only 0.23 residual standard deviations above the
  # trade-off line -- almost exactly where its coupling predicts. Being in a
  # rare corner of the map and being an outlier against the trade-off are
  # different things, which is precisely why the residual is a supplement to
  # the computed hierarchy and never a substitute for it.
  expect_equal(r$direction, "above")
  expect_equal(r$residual_sd, 0.226, tolerance = 1e-2)
  expect_lt(abs(r$residual_sd), 1)
  expect_match(r$caveat, "never instead of it")
})

test_that("everything returns NULL rather than failing on an undiagnosable matrix", {
  d <- spectral_diagnostics(matrix(2, 4, 4))   # uniform totals
  expect_true(is.na(d$mu_max))

  expect_null(structural_type(d))
  expect_null(type_stability(d))
  expect_null(tradeoff_residual(d))
  expect_null(structural_type(NULL))
})
