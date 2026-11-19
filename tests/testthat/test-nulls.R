# The two matched nulls: the full shuffle and the shuffle within rows.
#
# What this file catches: a null changing its distribution, the full shuffle
# drifting from what 0.6.0 drew, and the within-row shuffle losing a row total.
#
# The literals come from data-raw/reference_values.R, which runs the released
# 0.6.0 and the paper's R/baseline/spectral_baseline.R in a separate R process.
# The shuffles are sequences of integers and are the same on every platform.
# The diagnostics of each draw pass through an eigendecomposition, so the
# literals are compared at 1e-10: the spectra are well separated and another
# BLAS moves them in the last few digits only.

fefo <- worked_matrices$fefo_stock_control
tol <- 1e-10

test_that("a seeded full ensemble equals the 0.6.0 output", {
  # Released 0.6.0, surrogate_ensemble(fefo, B = 8, seed = 99).
  old <- list(
    mu_max = c(0.74839214937075504, 0.74985015876907657, 0.79282714142465027,
               0.75328060128263152, 0.77913526498664087, 0.77827825715245869,
               0.77844311745817418, 0.78227827449924447),
    lambda_max = c(2.9744387843984375, 2.9976039763977278, 3.8268871071076891,
                   3.0531875693550896, 3.5276580706264169, 3.5101575838126644,
                   3.5135135885982645, 3.5930188992394796),
    dominance = c(0.038398865471679741, 0.043744780428341359,
                  0.043600261394041377, 0.044635011092086938,
                  0.034156367968566873, 0.043318552800885188,
                  0.040832186591186287, 0.054830588148572112),
    hierarchy_sd = c(0.043109852051979668, 0.066590170053753867,
                     0.046623026034731029, 0.05561537241528345,
                     0.04754012568224375, 0.030973567164512824,
                     0.027708506278825852, 0.0460106334587434),
    hierarchy_pr = c(0.9131025446836627, 0.84923026358924036,
                     0.91080650064689717, 0.88583161233223739,
                     0.90832317393027417, 0.95671886823733854,
                     0.96286459367142052, 0.90649930787324251),
    hierarchy_gini = c(0.083348102004971988, 0.13107015221319376,
                       0.092055972814274453, 0.10894071690129079,
                       0.093954171241266218, 0.060726447647135723,
                       0.05347186594320627, 0.090654380212877861))

  ens <- surrogate_ensemble(fefo, B = 8, seed = 99)
  # The six columns of 0.6.0 come first and in the same order.
  expect_identical(names(ens)[1:6], names(old))
  for (q in names(old)) {
    expect_equal(ens[[q]], old[[q]], tolerance = tol, info = q)
  }
  expect_equal(attr(ens, "requested"), 8)
  expect_equal(attr(ens, "attempts"), 8)
  expect_true(attr(ens, "complete"))

  # Naming the default scheme changes nothing.
  expect_identical(ens, surrogate_ensemble(fefo, B = 8, seed = 99, scheme = "full"))
})

test_that("surrogate_position under reading sd and the full scheme equals 0.6.0", {
  # Released 0.6.0, surrogate_position(fefo, B = 25, seed = 99).
  sp <- surrogate_position(fefo, B = 25, seed = 99, reading = "sd")
  m <- sp$metrics[1:3, ]

  expect_identical(m$metric, c("mu_max", "hierarchy_sd", "dominance"))
  expect_identical(names(sp$metrics)[1:7],
                   c("metric", "observed", "median", "min", "max", "share_ge",
                     "outside"))
  expect_equal(m$observed, c(0.82783388004453395, 0.037756755992266425,
                             0.027887448063316951), tolerance = tol)
  expect_equal(m$median, c(0.76957297206838315, 0.04754012568224375,
                           0.043318552800885188), tolerance = tol)
  expect_equal(m$min, c(0.67773517575812947, 0.027708506278825852,
                        0.02913993235056898), tolerance = tol)
  expect_equal(m$max, c(0.80172132983327349, 0.066842813584860869,
                        0.059258710590174304), tolerance = tol)
  expect_equal(m$share_ge, c(0, 0.84, 1))
  expect_identical(m$outside, c(TRUE, FALSE, TRUE))

  expect_equal(sp$type, "diffuse-amplified")
  expect_equal(sp$type_share, 1)
  expect_equal(sp$B, 25)
  expect_equal(sp$attempts, 25)
  expect_equal(sp$hierarchy_reading, "hierarchy_sd")
  expect_equal(sp$scheme, "full")
})

test_that("both schemes reproduce the paper's nulls from the same seed", {
  # null_distribution() of the paper's code on the FEFO matrix, B = 6, after
  # set.seed(20261119). The package ports the two shuffles and draws them in
  # the same order, so the same seed gives the same six matrices.
  paper <- list(
    full = list(
      mu_max = c(0.80231020478382697, 0.73166729779311113, 0.70588871080763194,
                 0.65733363564876102, 0.75005687472288352, 0.67938708422224436),
      hierarchy_h = c(0.017880365133307929, 0.057062451437276951,
                      0.022916609732781534, 0.061452393326753008,
                      0.035949252238518747, 0.045748503495336323),
      dominance = c(0.040406468423167238, 0.046819979791533364,
                    0.043999934668232064, 0.061115362384650006,
                    0.047845019853661741, 0.056264294141262453)),
    rows = list(
      mu_max = c(0.69206490220195183, 0.78245353728155587, 0.74307617968555428,
                 0.76551921642819443, 0.68657228641248447, 0.72674770978136094),
      hierarchy_h = c(0.018424793221399782, 0.020724721752391928,
                      0.016808003344192524, 0.017688927522134024,
                      0.019599329386245645, 0.019947708274903886),
      dominance = c(0.050273439462165663, 0.036498837607668061,
                    0.047870893126695069, 0.046702301785246293,
                    0.06130888016859834, 0.048001146163035517)))

  for (scheme in names(paper)) {
    ens <- surrogate_ensemble(fefo, B = 6, seed = 20261119, scheme = scheme)
    expect_equal(attr(ens, "scheme"), scheme)
    for (q in names(paper[[scheme]])) {
      expect_equal(ens[[q]], paper[[scheme]][[q]], tolerance = tol,
                   info = paste(scheme, q))
    }
  }
})

test_that("the ensemble mean, standard deviation and z match the paper's", {
  # compare_to_null() of the paper's code on the six draws of the rows null
  # above: null_mean, null_sd and z.
  sp <- surrogate_position(fefo, B = 6, seed = 20261119, scheme = "rows")
  m <- sp$metrics
  expect_identical(m$metric, c("mu_max", "hierarchy_sd", "dominance", "hierarchy_h"))
  expect_identical(names(m), c("metric", "observed", "median", "min", "max",
                               "share_ge", "outside", "mean", "sd", "z"))

  ref <- data.frame(
    metric = c("mu_max", "dominance", "hierarchy_sd", "hierarchy_h"),
    mean = c(0.73273897196518367, 0.048442583052234818, 0.039624422909667729,
             0.018865580583544633),
    sd = c(0.038664139260265538, 0.0079384755049225798, 0.0015633903048896183,
           0.00148177080231917),
    z = c(2.4595118344475626, -2.5893050846062455, -1.1946261349837202,
          -1.1868968561654167))
  for (i in seq_len(nrow(ref))) {
    row <- m[m$metric == ref$metric[i], ]
    expect_equal(row$mean, ref$mean[i], tolerance = tol, info = ref$metric[i])
    # The standard deviation of six close numbers loses digits to
    # cancellation, and z divides by it.
    expect_equal(row$sd, ref$sd[i], tolerance = 1e-7, info = ref$metric[i])
    expect_equal(row$z,  ref$z[i],  tolerance = 1e-7, info = ref$metric[i])
  }
  expect_equal(sp$scheme, "rows")
  expect_equal(sp$hierarchy_reading, "hierarchy_h")
})

test_that("z is the observed value less the mean, over the standard deviation", {
  sp <- surrogate_position(fefo, B = 30, seed = 4)
  m <- sp$metrics
  expect_equal(m$z, (m$observed - m$mean) / m$sd)
  ens <- surrogate_ensemble(fefo, B = 30, seed = 4)
  expect_equal(m$mean, vapply(m$metric, function(q) mean(ens[[q]]), numeric(1)),
               ignore_attr = TRUE)
  expect_equal(m$sd, vapply(m$metric, function(q) stats::sd(ens[[q]]), numeric(1)),
               ignore_attr = TRUE)
})

test_that("every rows draw keeps each row total and each row's multiset", {
  shuffle <- spectralDEMATEL:::shuffle_matrix
  set.seed(31)
  for (A in list(fefo, random_admissible(6), random_admissible(11))) {
    attributes(A) <- list(dim = dim(A))
    for (draw in 1:40) {
      S <- shuffle(A, "rows")
      expect_identical(dim(S), dim(A))
      expect_identical(diag(S), diag(A))
      # Sums of the same numbers in another order, so equal up to rounding.
      expect_equal(rowSums(S), rowSums(A), tolerance = 1e-13)
      for (i in seq_len(nrow(A))) {
        expect_identical(sort(S[i, -i]), sort(A[i, -i]))
      }
    }
  }
})

test_that("the rows null moves the column totals, so it is a null", {
  # If the shuffle returned the matrix unchanged the invariants above would
  # hold trivially.
  shuffle <- spectralDEMATEL:::shuffle_matrix
  set.seed(32)
  A <- fefo
  moved <- vapply(1:20, function(i) {
    S <- shuffle(A, "rows")
    !isTRUE(all.equal(colSums(S), colSums(A)))
  }, logical(1))
  expect_true(all(moved))
})

test_that("every full draw keeps the multiset of off-diagonal entries", {
  shuffle <- spectralDEMATEL:::shuffle_matrix
  set.seed(33)
  A <- fefo
  attributes(A) <- list(dim = dim(A))
  off <- row(A) != col(A)
  changed <- 0
  for (draw in 1:20) {
    S <- shuffle(A, "full")
    expect_identical(diag(S), diag(A))
    expect_identical(sort(S[off]), sort(A[off]))
    if (!isTRUE(all.equal(rowSums(S), rowSums(A)))) changed <- changed + 1
  }
  # The full shuffle does not keep the row totals.
  expect_gt(changed, 15)
})

test_that("the ensemble is the shuffle applied draw by draw", {
  # Ties the invariants above, asserted on the internal shuffle, to what the
  # exported function returns.
  shuffle <- spectralDEMATEL:::shuffle_matrix
  for (scheme in c("full", "rows")) {
    ens <- surrogate_ensemble(fefo, B = 5, seed = 77, scheme = scheme)
    set.seed(77)
    by_hand <- vapply(1:5, function(i) {
      spectral_diagnostics(shuffle(fefo, scheme), checks = FALSE)$mu_max
    }, numeric(1))
    expect_equal(ens$mu_max, by_hand, tolerance = 1e-14, info = scheme)
  }
})

test_that("hierarchy_h in the ensemble follows its identity with sd", {
  for (scheme in c("full", "rows")) {
    ens <- surrogate_ensemble(fefo, B = 20, seed = 6, scheme = scheme)
    expect_identical(names(ens), c("mu_max", "lambda_max", "dominance",
                                   "hierarchy_sd", "hierarchy_pr",
                                   "hierarchy_gini", "hierarchy_h"))
    expect_equal(ens$hierarchy_h, (nrow(fefo) - 1) * ens$hierarchy_sd^2,
                 tolerance = 1e-12)
  }
})

test_that("a rows null with nothing to shuffle has a zero spread and z of NA", {
  # Every row is constant off the diagonal, so a shuffle within rows returns
  # the matrix itself. The spread is zero and z is undefined. The paper's
  # compare_to_null() returns NA there and so does this.
  A <- outer(c(1, 2, 3, 4), rep(1, 4))
  diag(A) <- 0
  expect_true(is_irreducible(A))

  sp <- surrogate_position(A, B = 10, seed = 1, scheme = "rows")
  expect_equal(sp$metrics$sd, rep(0, 4), tolerance = 1e-12)
  expect_true(all(is.na(sp$metrics$z)))
  expect_equal(sp$metrics$mean, sp$metrics$observed, tolerance = 1e-12)
})

test_that("rejection under the rows scheme is bounded like the full scheme", {
  # A weighted directed cycle. Each row holds one rating, and a shuffle within
  # rows keeps the graph strongly connected only when the five ratings form a
  # single cycle again: 24 of the 1024 arrangements.
  n <- 5
  cyc <- matrix(0, n, n)
  for (i in 1:n) cyc[i, if (i == n) 1 else i + 1] <- i
  expect_true(is_irreducible(cyc))

  ens <- surrogate_ensemble(cyc, B = 100, seed = 1, max_attempts = 300,
                            scheme = "rows")
  expect_s3_class(ens, "data.frame")
  expect_lt(nrow(ens), 100)
  expect_false(attr(ens, "complete"))
  expect_equal(attr(ens, "attempts"), 300)
  expect_equal(attr(ens, "scheme"), "rows")

  sp <- surrogate_position(cyc, B = 100, seed = 1, max_attempts = 300,
                           scheme = "rows")
  expect_false(sp$complete)
  expect_equal(sp$B, nrow(ens))
})

test_that("a two-factor system has a rows null, though nothing can move", {
  ens <- surrogate_ensemble(two_by_two, B = 4, seed = 1, scheme = "rows")
  expect_equal(nrow(ens), 4)
  expect_equal(ens$mu_max,
               rep(spectral_diagnostics(two_by_two, checks = FALSE)$mu_max, 4),
               tolerance = 1e-12)
})

test_that("the type share is counted on the reading that was asked for", {
  sp_h  <- surrogate_position(fefo, B = 20, seed = 9)
  sp_sd <- surrogate_position(fefo, B = 20, seed = 9, reading = "sd")
  expect_equal(sp_h$hierarchy_reading, "hierarchy_h")
  expect_equal(sp_sd$hierarchy_reading, "hierarchy_sd")
  # Same draws, so the metric rows agree whatever the reading.
  expect_equal(sp_h$metrics, sp_sd$metrics)
  expect_equal(sp_h$type, structural_type(spectral_diagnostics(fefo))$type)
})
