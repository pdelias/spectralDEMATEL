# Hierarchy as H = 1 - Neff / n, and the two readings that classify.
#
# What this file catches: H drifting from its identity with sd(u), H running
# the wrong way, and reading = "sd" no longer reproducing 0.6.0.
#
# The independent routes are the paper's own code run in a separate R process
# (data-raw/reference_values.R), the identity H = (n - 1) sd(u)^2, the replica
# of Appendix A of the source paper, and the decision model of its Appendix B
# evaluated numerically.

fefo <- worked_matrices$fefo_stock_control
resl <- worked_matrices$resilience_capabilities

test_that("H of the two worked matrices matches the paper's code", {
  # hierarchy_H from diagnostics() in the paper's R/baseline/spectral_baseline.R,
  # printed to 17 significant digits by data-raw/reference_values.R.
  #
  # Tolerance. Both routes read the same transcribed matrix, so the only
  # difference is rounding. H is one minus a ratio close to one, so its
  # absolute error is a small multiple of machine epsilon times the eigenvector
  # condition. Both matrices have a well separated dominant eigenvalue
  # (dominance 0.028 and 0.199) and a condition number of 1.01 and 2.17. The
  # measured differences were 4.4e-16 and 2.2e-16. An absolute bound of 1e-12
  # leaves more than three orders of magnitude for another BLAS.
  expect_lt(abs(spectral_diagnostics(fefo)$hierarchy_h - 0.017106871476714303), 1e-12)
  expect_lt(abs(spectral_diagnostics(resl)$hierarchy_h - 0.3056477794475575), 1e-12)
})

test_that("H matches the values the manuscript prints", {
  # Keys quad_worked_hier_lo and quad_worked_hier_hi of the paper's results,
  # rounded there to three decimals. The pipeline zeroes the diagonal and this
  # package keeps it. For these two matrices H is the same either way: the FEFO
  # diagonal is constant and the resilience diagonal is zero.
  expect_equal(round(spectral_diagnostics(fefo)$hierarchy_h, 3), 0.017)
  expect_equal(round(spectral_diagnostics(resl)$hierarchy_h, 3), 0.306)
})

test_that("H equals (n - 1) sd(u)^2 on every fixture", {
  set.seed(11)
  fixtures <- list(fefo = fefo, resl = resl, cycle3 = cycle3,
                   defective = defective, two_by_two = two_by_two,
                   rank1 = rank1,
                   random6 = random_admissible(6),
                   random12 = random_admissible(12),
                   random25 = random_admissible(25))
  for (nm in names(fixtures)) {
    d <- spectral_diagnostics(fixtures[[nm]], checks = FALSE)
    expect_false(is.null(d), info = nm)
    # Both sides are computed from the same vector, so this is rounding only.
    expect_lt(abs(d$hierarchy_h - (d$n - 1) * d$hierarchy_sd^2), 1e-12)
    # The definition, written out from the entry profile.
    u <- d$entry_points
    expect_lt(abs(d$hierarchy_h - (1 - sum(u)^2 / (d$n * sum(u^2)))), 1e-12)
  }
})

test_that("H is high when influence enters at a few factors", {
  # A hub that dispatches to every factor and receives a little from each,
  # against a matrix with nearly equal ratings. The direction is the point.
  n <- 8
  hub <- matrix(0, n, n)
  hub[1, -1] <- 4
  hub[-1, 1] <- 0.2
  even <- matrix(2, n, n)
  diag(even) <- 0
  even[1, 2] <- 2.1

  h_hub  <- spectral_diagnostics(hub)$hierarchy_h
  h_even <- spectral_diagnostics(even)$hierarchy_h
  expect_gt(h_hub, 0.3)
  expect_lt(h_even, 0.001)

  # It runs with the standard deviation and the Gini coefficient and against
  # the participation ratio.
  d_hub  <- spectral_diagnostics(hub)
  d_even <- spectral_diagnostics(even)
  expect_gt(d_hub$hierarchy_sd,   d_even$hierarchy_sd)
  expect_gt(d_hub$hierarchy_gini, d_even$hierarchy_gini)
  expect_lt(d_hub$hierarchy_pr,   d_even$hierarchy_pr)
})

test_that("H stays between zero and 1 - 1/n", {
  set.seed(3)
  for (n in c(2, 5, 9, 20)) {
    d <- spectral_diagnostics(random_admissible(n), checks = FALSE)
    expect_gte(d$hierarchy_h, 0)
    expect_lt(d$hierarchy_h, 1 - 1 / n)
  }
})

test_that("H and the Gini coefficient do not change when a vector is replicated", {
  # Ported from section 3 of the paper's R/baseline/test_spectral_baseline.R.
  # The paper's evenness(g) is sum(g)^2 / (length(g) * sum(g^2)), and H is one
  # minus it.
  evenness <- function(g) sum(g)^2 / (length(g) * sum(g^2))
  set.seed(20261024)
  for (trial in 1:200) {
    n <- sample(4:30, 1)
    v <- stats::rexp(n)^stats::runif(1, 0.2, 3)
    v <- v / sqrt(sum(v^2))
    expect_lt(abs((n - 1) * stats::sd(v)^2 - (1 - evenness(v))), 1e-12)
    k <- sample(2:4, 1)
    v_rep <- rep(v, each = k) / sqrt(k)
    expect_lt(abs(evenness(v_rep) - evenness(v)), 1e-12)
    expect_lt(abs(spectralDEMATEL:::gini(v_rep) - spectralDEMATEL:::gini(v)), 1e-12)
    expect_lt(stats::sd(v_rep), stats::sd(v))
  }
})

test_that("a replica keeps coupling, mode dominance and H", {
  # Appendix A of the source paper. The replica of a model with direct matrix A
  # is the Kronecker product of A with the k by k matrix of ones: each factor
  # becomes k copies with the same ratings. Its spectrum is that of the model
  # with n(k - 1) zeros added, so coupling, mode dominance and H stay, while
  # sd(u)^2 falls by the factor (n - 1) / (kn - 1).
  #
  # Tolerance. The added eigenvalue is zero with multiplicity n(k - 1), far
  # from the dominant one, so the dominant pair is as well conditioned as in
  # the model. 1e-9 covers the larger matrix.
  set.seed(8)
  models <- list(fefo = fefo, random7 = random_admissible(7),
                 random10 = random_admissible(10))
  for (nm in names(models)) {
    A <- models[[nm]]
    attributes(A) <- list(dim = dim(A))
    d <- spectral_diagnostics(A, checks = FALSE)
    for (k in 2:3) {
      R <- kronecker(A, matrix(1, k, k))
      r <- spectral_diagnostics(R, checks = FALSE)
      info <- paste(nm, "k =", k)
      expect_equal(r$n, k * d$n, info = info)
      expect_lt(abs(r$mu_max - d$mu_max), 1e-9)
      expect_lt(abs(r$dominance - d$dominance), 1e-9)
      expect_lt(abs(r$hierarchy_h - d$hierarchy_h), 1e-9)
      expect_lt(abs(r$hierarchy_sd^2 -
                      d$hierarchy_sd^2 * (d$n - 1) / (k * d$n - 1)), 1e-9)
      # The Gini coefficient is replication invariant too.
      expect_lt(abs(r$hierarchy_gini - d$hierarchy_gini), 1e-9)
    }
  }
})

test_that("an inadmissible matrix has hierarchy_h in its NA shape", {
  bad <- spectral_diagnostics(matrix(2, 4, 4))
  expect_true("hierarchy_h" %in% names(bad))
  expect_true(is.na(bad$hierarchy_h))
})

test_that("a published total-relation matrix gives the same H", {
  m <- dematel(fefo)
  from_T <- spectral_diagnostics(m$T, type = "T", checks = FALSE)
  from_A <- spectral_diagnostics(fefo, checks = FALSE)
  expect_lt(abs(from_T$hierarchy_h - from_A$hierarchy_h), 1e-10)
})

# ---- the readings -----------------------------------------------------------

# A diagnostics list with a chosen coupling, sd and size, and the H that the
# identity gives for them.
mk <- function(mu, sd, n) {
  list(n = n, mu_max = mu, hierarchy_sd = sd, hierarchy_h = (n - 1) * sd^2)
}

test_that("the default reading is H, cut at 0.10", {
  s <- structural_type(spectral_diagnostics(fefo))
  expect_equal(s$hierarchy_reading, "hierarchy_h")
  expect_match(s$hierarchy_direction, "high")
  expect_equal(s$cuts$hierarchy, 0.10)
  expect_true(s$cuts$hierarchy_is_default)
  expect_equal(s$hierarchy_margin,
               spectral_diagnostics(fefo)$hierarchy_h - 0.10)
  expect_equal(s$type, "diffuse-amplified")

  r <- structural_type(spectral_diagnostics(resl))
  expect_equal(r$type, "hierarchical-dampened")
  expect_gt(r$hierarchy_margin, 0)
})

test_that("reading = \"sd\" reproduces the 0.6.0 classification", {
  s <- structural_type(spectral_diagnostics(fefo), reading = "sd")
  expect_equal(s$hierarchy_reading, "hierarchy_sd")
  expect_equal(s$type, "diffuse-amplified")
  expect_equal(s$corpus_share, 81 / 117)
  expect_equal(s$corpus_multiplier, 5.14)
  expect_equal(s$corpus_n, 117L)
  expect_equal(s$hierarchy_margin, spectral_diagnostics(fefo)$hierarchy_sd - 0.10)
  # 0.6.0 scaled the margins by the interquartile ranges of the 117 systems.
  expect_equal(s$nearest_margin_scaled,
               min(abs(s$coupling_margin) / (0.8555 - 0.5824),
                   abs(s$hierarchy_margin) / (0.0985 - 0.0279)))
})

test_that("the two readings place the cut differently, except at eleven factors", {
  # For u of unit length H = (n - 1) sd^2, so sd = 0.10 is H = 0.10 only when
  # n = 11. A large model with a modest sd is hierarchical on H and diffuse on
  # sd, and a small model the other way round.
  large <- mk(0.8, 0.07, 30)        # H = 0.1421
  expect_equal(structural_type(large)$type, "hierarchical-amplified")
  expect_equal(structural_type(large, reading = "sd")$type, "diffuse-amplified")

  small <- mk(0.8, 0.12, 5)         # H = 0.0576
  expect_equal(structural_type(small)$type, "diffuse-amplified")
  expect_equal(structural_type(small, reading = "sd")$type, "hierarchical-amplified")

  for (sd in c(0.05, 0.099, 0.101, 0.2)) {
    eleven <- mk(0.8, sd, 11)
    expect_equal(structural_type(eleven)$type,
                 structural_type(eleven, reading = "sd")$type)
  }
})

test_that("a reading whose value is missing returns NULL", {
  # A diagnostics list saved by 0.6.0 has no hierarchy_h. It can still be
  # classified on sd. It is not silently reclassified on a derived H.
  old <- list(n = 13, mu_max = 0.8, hierarchy_sd = 0.05)
  expect_null(structural_type(old))
  expect_equal(structural_type(old, reading = "sd")$type, "diffuse-amplified")
  expect_null(type_stability(old))
  expect_null(tradeoff_residual(old))
})

test_that("the corpus behind each reading is the one that reading was cut on", {
  s_h  <- structural_type(spectral_diagnostics(fefo))
  s_sd <- structural_type(spectral_diagnostics(fefo), reading = "sd")
  cp <- spectralDEMATEL:::CORPUS
  expect_equal(s_h$corpus_n,  cp$H$n)
  expect_equal(s_sd$corpus_n, cp$sd$n)
  expect_equal(s_h$corpus_share,
               unname(cp$H$count["diffuse-amplified"]) / cp$H$n)
  expect_equal(s_h$corpus_multiplier,
               unname(cp$H$multiplier["diffuse-amplified"]))
})

test_that("type_stability uses the grid of its reading", {
  d <- spectral_diagnostics(fefo)
  st <- type_stability(d)
  expect_equal(st$by_cut$cut, c(0.05, 0.075, 0.10, 0.125, 0.20))
  expect_equal(st$hierarchy_reading, "hierarchy_h")
  expect_true(st$stable)

  st_sd <- type_stability(d, reading = "sd")
  expect_equal(st_sd$by_cut$cut, c(0.075, 0.09, 0.10, 0.11))
  expect_equal(st_sd$hierarchy_reading, "hierarchy_sd")

  # A model inside the band flips, and says between which two cuts.
  inside <- mk(0.8, sqrt(0.11 / 12), 13)     # H = 0.11
  st2 <- type_stability(inside)
  expect_false(st2$stable)
  expect_equal(st2$flips_between, c(0.10, 0.125))
})

test_that("the trade-off residual is read on the corpus of its reading", {
  d <- spectral_diagnostics(resl)
  cp <- spectralDEMATEL:::CORPUS

  r <- tradeoff_residual(d)
  expect_equal(r$hierarchy_reading, "hierarchy_h")
  expect_equal(r$expected, cp$H$tradeoff$intercept + cp$H$tradeoff$slope * d$mu_max)
  expect_equal(r$residual, d$hierarchy_h - r$expected)
  expect_equal(r$residual_sd, r$residual / cp$H$tradeoff$resid_sd)
  expect_equal(r$fit, cp$H$tradeoff$method)
  expect_match(r$caveat, "9 of the 81")

  r_sd <- tradeoff_residual(d, reading = "sd")
  expect_equal(r_sd$hierarchy_reading, "hierarchy_sd")
  expect_equal(r_sd$residual_sd, 0.226, tolerance = 1e-2)
  expect_match(r_sd$caveat, "19 of the 117")
})

# ---- the effort shares --------------------------------------------------------

test_that("the effort shares are H and the square of coupling", {
  d <- spectral_diagnostics(fefo)
  e <- structural_type(d)$effort
  expect_named(e, c("uniform_forgoes", "direct_appraisal_forgoes", "assumptions"))
  expect_equal(e$uniform_forgoes, d$hierarchy_h)
  expect_equal(e$direct_appraisal_forgoes, d$mu_max^2)
  expect_match(e$assumptions, "equal costs")
  expect_match(e$assumptions, "dominant mode")

  # One quarter at the coupling cut.
  expect_equal(structural_type(mk(0.5, 0.05, 10))$effort$direct_appraisal_forgoes, 0.25)

  # Under the sd reading the shares are not reported.
  e_sd <- structural_type(d, reading = "sd")$effort
  expect_named(e_sd, c("uniform_forgoes", "direct_appraisal_forgoes", "assumptions"))
  expect_true(is.na(e_sd$uniform_forgoes))
  expect_true(is.na(e_sd$direct_appraisal_forgoes))
})

test_that("the effort shares agree with the decision model evaluated numerically", {
  # Appendix B of the source paper. Effort x has net value
  #   W(x) = x'g - x'Cx / 2,  g = (I + T) w,  C = diag(costs),
  # and the dominant-mode objective sets w = u. The best allocation is
  # C^-1 g. Uniform effort is t 1 with the best t. Direct appraisal acts as if
  # T were zero and chooses C^-1 u. Nothing below uses mu_max or H: the shares
  # are computed from T, u and W alone.
  set.seed(5)
  for (A in list(fefo, random_admissible(7), random_admissible(15))) {
    attributes(A) <- list(dim = dim(A))
    d <- spectral_diagnostics(A, checks = FALSE)
    Tm <- dematel(A)$T
    n <- d$n
    u <- d$entry_points
    g <- as.vector((diag(n) + Tm) %*% u)
    e <- structural_type(d)$effort

    W <- function(x, cost) sum(x * g) - sum(cost * x^2) / 2

    # Uniform effort, equal costs.
    cost <- rep(1, n)
    best <- W(g / cost, cost)
    t_star <- sum(g) / sum(cost)
    expect_lt(abs((1 - W(rep(t_star, n), cost) / best) - e$uniform_forgoes), 1e-10)

    # Direct appraisal, whatever the costs.
    for (trial in 1:3) {
      cost <- stats::runif(n, 0.5, 3)
      best <- W(g / cost, cost)
      expect_lt(abs((1 - W(u / cost, cost) / best) - e$direct_appraisal_forgoes),
                1e-10)
    }
  }
})

test_that("robustness functions report the reading they classified on", {
  ms <- measurement_stability(fefo, B = 10, seed = 2)
  expect_equal(ms$hierarchy_reading, "hierarchy_h")
  ms_sd <- measurement_stability(fefo, B = 10, seed = 2, reading = "sd")
  expect_equal(ms_sd$hierarchy_reading, "hierarchy_sd")

  # Released 0.6.0, measurement_stability(resl, B = 40, seed = 5).
  old <- measurement_stability(resl, B = 40, seed = 5, reading = "sd")
  expect_equal(old$observed_type, "hierarchical-dampened")
  expect_equal(old$share_same, 1)
  expect_equal(old$B_admissible, 40)
})

# ---- the corpus constants of the H reading -----------------------------------

test_that("the H corpus constants are internally consistent", {
  cp <- spectralDEMATEL:::CORPUS$H

  expect_equal(cp$n, 81L)
  expect_equal(sum(cp$count), cp$n)
  expect_setequal(names(cp$count), names(cp$multiplier))
  expect_setequal(names(cp$count),
                  names(spectralDEMATEL:::INTERVENTION_LOGIC$H))

  # The guide to the panel method and the manuscript both give 57
  # diffuse-amplified models.
  expect_equal(unname(cp$count["diffuse-amplified"]), 57L)
  # Shares the manuscript's Table of types prints, to the nearest per cent.
  share <- round(100 * cp$count / cp$n)
  expect_equal(unname(share["diffuse-amplified"]),      70)
  expect_equal(unname(share["hierarchical-dampened"]),  16)
  expect_equal(unname(share["hierarchical-amplified"]), 10)
  expect_equal(unname(share["diffuse-dampened"]),        4)

  expect_named(cp$quantiles, c("mu_max", "hierarchy_h", "dominance"))
  for (q in cp$quantiles) {
    expect_identical(q, sort(q))
    expect_length(q, 5)
  }
  expect_length(cp$probs, 5)

  # The upper quartile of H lies just above the cut, as the manuscript says
  # (key quad_hier_q75 = 0.11).
  expect_equal(round(cp$quantiles$hierarchy_h[4], 2), 0.11)

  expect_lt(cp$tradeoff$slope, 0)
  expect_gt(cp$tradeoff$resid_sd, 0)
  expect_equal(round(cp$tradeoff$r2, 2), 0.75)
  expect_equal(cp$tradeoff$spearman, -0.879)
  expect_true(nzchar(cp$tradeoff$method))
  expect_true(nzchar(cp$source$generated_at))
})

test_that("the H trade-off line reproduces the counts the manuscript reports", {
  # An independent route to the intercept and the slope. Used alone, the line
  # assigns a type from coupling. The manuscript reports that this is wrong for
  # every diffuse-dampened model: at a coupling below one half the line
  # predicts a hierarchy above the cut. So the line must cross H = 0.10 at a
  # coupling above 0.5.
  tr <- spectralDEMATEL:::CORPUS$H$tradeoff
  crossing <- (0.10 - tr$intercept) / tr$slope
  expect_gt(crossing, 0.5)
  expect_lt(crossing, 1)
})

test_that("the sd corpus constants are those of 0.6.0", {
  cp <- spectralDEMATEL:::CORPUS$sd
  expect_equal(cp$n, 117L)
  expect_equal(cp$tradeoff$intercept, 0.211766)
  expect_equal(cp$tradeoff$slope, -0.204723)
  expect_equal(cp$tradeoff$resid_sd, 0.034985)
  expect_equal(unname(cp$count), c(81L, 16L, 12L, 8L))
  expect_equal(unname(cp$multiplier), c(5.14, 1.50, 3.06, 1.83))
  expect_equal(cp$quantiles$mu_max, c(0.2432, 0.5824, 0.7541, 0.8555, 0.9187))
  expect_equal(cp$quantiles$hierarchy_sd, c(0.0140, 0.0279, 0.0549, 0.0985, 0.1686))
  expect_equal(cp$quantiles$dominance, c(0.0121, 0.0253, 0.0457, 0.1032, 0.3389))
})
