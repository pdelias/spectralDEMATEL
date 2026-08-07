# The assumption checks.
#
# Two things matter more than any individual verdict. The table must be
# rectangular whatever is passed in, because a batch job rbinds a hundred of
# them. And `skipped` must never be confusable with `pass`, because an interface
# that renders them the same way tells a user their matrix is fine when nothing
# was actually tested.

test_that("the table is rectangular and complete for every input", {
  inputs <- list(
    good         = worked_matrices$fefo_stock_control,
    disconnected = worked_matrices$resilience_capabilities,
    not_square   = matrix(1:6, 2, 3),
    not_a_matrix = as.data.frame(matrix(1, 2, 2)),
    has_na       = matrix(c(0, NA, 1, 0), 2, 2),
    negative     = matrix(c(0, -1, 1, 0), 2, 2),
    all_zero     = matrix(0, 4, 4),
    uniform      = matrix(2, 4, 4),
    character    = matrix(letters[1:4], 2, 2)
  )

  for (nm in names(inputs)) {
    ck <- assumption_checks(inputs[[nm]])
    expect_s3_class(ck, "spectral_checks")
    expect_s3_class(ck, "data.frame")
    expect_named(ck, c("check", "verdict", "reason", "value", "factors"))
    expect_identical(ck$check, spectralDEMATEL:::CHECK_IDS, info = nm)
    expect_true(all(ck$verdict %in% c("pass", "warn", "fail", "skipped")),
                info = nm)
    expect_true(all(nzchar(ck$reason)), info = nm)
    expect_type(ck$factors, "list")
  }
})

test_that("a hundred matrices stack into one table", {
  # The batch-mode requirement, asserted rather than assumed.
  set.seed(8)
  many <- lapply(1:20, function(i) {
    if (i %% 5 == 0) matrix(2, 4, 4) else random_admissible(sample(4:9, 1))
  })
  stacked <- do.call(rbind, lapply(many, assumption_checks))
  expect_equal(nrow(stacked), 20 * length(spectralDEMATEL:::CHECK_IDS))
})

test_that("a clean matrix passes everything it should", {
  ck <- assumption_checks(worked_matrices$fefo_stock_control)
  verdict <- function(id) ck$verdict[ck$check == id]

  expect_equal(verdict("input_is_matrix"),     "pass")
  expect_equal(verdict("square"),              "pass")
  expect_equal(verdict("finite"),              "pass")
  expect_equal(verdict("nonnegative"),         "pass")
  expect_equal(verdict("has_influence"),       "pass")
  expect_equal(verdict("strong_connectivity"), "pass")
  expect_equal(verdict("totals_vary"),         "pass")
  expect_equal(verdict("connectivity_margin"), "pass")
  expect_equal(verdict("coupling_margin"),     "pass")   # 0.828, below 0.95

  # But its diagonal is not zero, which is admissible and worth saying.
  expect_equal(verdict("zero_diagonal"), "warn")
  expect_equal(ck$value[ck$check == "zero_diagonal"], 13)
})

test_that("the disconnected worked matrix names the factors at fault", {
  # This published matrix strands four of its fifteen factors, and the chain of
  # cause is worth spelling out because it is what a user needs to see:
  #
  #   14  row 14 is all zeros, so it dispatches nothing
  #    7  column 7 is all zeros, so it receives nothing
  #    8  its only source is factor 7, which is itself unreachable
  #    4  its only sources are 7 and 8, likewise
  #
  # Naming these four is the point. A user is told where to look.
  A <- worked_matrices$resilience_capabilities
  expect_identical(which(rowSums(A) == 0), 14L)
  expect_identical(which(colSums(A) == 0), 7L)
  expect_identical(which(A[, 8] > 0), 7L)
  expect_identical(which(A[, 4] > 0), c(7L, 8L))

  ck <- assumption_checks(A)
  expect_equal(ck$verdict[ck$check == "strong_connectivity"], "fail")
  stranded <- ck$factors[ck$check == "strong_connectivity"][[1]]

  # Exactly those four. The naive definition -- everything that cannot reach or
  # be reached by all others -- names all fifteen, because once one factor is
  # stranded nobody can reach it either. That is useless to a user looking for
  # the row they mistyped.
  expect_identical(stranded, c(4L, 7L, 8L, 14L))
  expect_equal(ck$value[ck$check == "strong_connectivity"], 4)
  expect_match(ck$reason[ck$check == "strong_connectivity"], "cut off")

  # And the margin check cannot be evaluated on an already-broken graph.
  expect_equal(ck$verdict[ck$check == "connectivity_margin"], "skipped")
})

test_that("the reachability closure has exactly one implementation", {
  # is_irreducible() and the strong_connectivity check must never disagree,
  # which is guaranteed by their sharing reachability(). This is the test that
  # would fail if someone re-implemented one of them.
  set.seed(21)
  cases <- c(list(cycle3, defective, block_diagonal,
                  worked_matrices$fefo_stock_control,
                  worked_matrices$resilience_capabilities),
             lapply(4:9, random_admissible))

  for (A in cases) {
    ck <- assumption_checks(A)
    from_check <- ck$verdict[ck$check == "strong_connectivity"] == "pass"
    expect_identical(from_check, is_irreducible(A))

    # A connected graph strands nobody, and a broken one strands somebody.
    stranded <- ck$factors[ck$check == "strong_connectivity"][[1]]
    expect_identical(length(stranded) == 0L, is_irreducible(A))
  }
})

test_that("two separate components name one side, not both", {
  # block_diagonal is two 2-cycles that never meet. Naming all four would be
  # true and unhelpful; naming the half outside the main body is actionable.
  ck <- assumption_checks(block_diagonal)
  stranded <- ck$factors[ck$check == "strong_connectivity"][[1]]
  expect_length(stranded, 2)
  expect_true(setequal(stranded, c(3L, 4L)) || setequal(stranded, c(1L, 2L)))
})

test_that("a skipped check carries no value and no factors", {
  # Found in the application: connectivity_margin was reporting the factors
  # hanging on a single link even when the graph was already disconnected and
  # the check had never run. Names beside a verdict that was not reached read
  # as evidence for it.
  inputs <- list(worked_matrices$resilience_capabilities,
                 block_diagonal,
                 matrix(1:6, 2, 3),
                 matrix(2, 4, 4),
                 matrix(c(0, NA, 1, 0), 2, 2))

  for (A in inputs) {
    ck <- assumption_checks(A)
    skipped <- ck[ck$verdict == "skipped", ]
    if (nrow(skipped) == 0L) next
    expect_true(all(is.na(skipped$value)))
    expect_true(all(lengths(skipped$factors) == 0L))
  }
})

test_that("skipped propagates from the first structural failure onward", {
  ck <- assumption_checks(matrix(1:6, 2, 3))   # not square

  expect_equal(ck$verdict[ck$check == "input_is_matrix"], "pass")
  expect_equal(ck$verdict[ck$check == "square"], "fail")
  # Everything after `square` could not be evaluated. None of it is a pass.
  later <- ck$verdict[match("finite", ck$check):nrow(ck)]
  expect_true(all(later == "skipped"))
  expect_false(any(later == "pass"))
})

test_that("a missing entry blocks the sign check rather than passing it", {
  # any(A < 0) is NA when A holds NA, so a naive implementation reports a pass.
  ck <- assumption_checks(matrix(c(0, NA, 1, 0), 2, 2))
  expect_equal(ck$verdict[ck$check == "finite"], "fail")
  expect_equal(ck$verdict[ck$check == "nonnegative"], "skipped")
})

test_that("near-violations warn where a binary verdict would pass", {
  # Coupling: in scope, but close enough to criticality that the multiplier is
  # very sensitive. Built by pushing one row total close to the maximum.
  near <- matrix(c(0, 100, 1,
                   100, 0, 1,
                   1,   1, 0), 3, 3, byrow = TRUE)
  storage.mode(near) <- "double"
  d <- spectral_diagnostics(near)
  skip_if(is.na(d$mu_max))
  if (d$mu_max >= 0.95) {
    expect_equal(d$checks$verdict[d$checks$check == "coupling_margin"], "warn")
  }

  # Connectivity margin: strongly connected, but one factor hangs on a single
  # link, so removing one rating would break it.
  fragile <- matrix(c(0, 1, 1,
                      1, 0, 1,
                      1, 0, 0), 3, 3, byrow = TRUE)
  storage.mode(fragile) <- "double"
  ck <- assumption_checks(fragile)
  expect_equal(ck$verdict[ck$check == "strong_connectivity"], "pass")
  expect_equal(ck$verdict[ck$check == "connectivity_margin"], "warn")
  expect_gt(length(ck$factors[ck$check == "connectivity_margin"][[1]]), 0)
})

test_that("conditioning warns for a long chain and passes for a hub", {
  # Measured before the threshold was chosen: a symmetric hub gives exactly 1,
  # random admissible matrices run to 2, and a directed chain reaches 30.
  chain <- matrix(0, 10, 10)
  for (i in 1:9) chain[i, i + 1] <- 4
  chain[10, 1] <- 0.01
  d <- spectral_diagnostics(chain)
  expect_gt(d$ev_condition, spectralDEMATEL:::CHECK_THRESHOLDS$ev_condition_warn)
  expect_equal(d$checks$verdict[d$checks$check == "sensitivity_conditioning"],
               "warn")

  hub <- matrix(0, 8, 8)
  hub[1, 2:8] <- 4
  hub[2:8, 1] <- 4
  d2 <- spectral_diagnostics(hub)
  expect_equal(d2$ev_condition, 1, tolerance = 1e-8)
  expect_equal(d2$checks$verdict[d2$checks$check == "sensitivity_conditioning"],
               "pass")
})

test_that("the normalisation checks do not apply to a total-relation matrix", {
  A <- worked_matrices$fefo_stock_control
  ck <- assumption_checks(dematel(A)$T, type = "T")

  expect_equal(ck$verdict[ck$check == "totals_vary"], "skipped")
  expect_match(ck$reason[ck$check == "totals_vary"], "not normalised")
  # The spectral checks still run, because the spectrum is recoverable.
  expect_true(ck$verdict[ck$check == "coupling_margin"] %in% c("pass", "warn"))
})

test_that("reasons are written for a user, not a developer", {
  # No function names, no R types, no variable names leaking into text a
  # person is meant to read.
  for (A in list(worked_matrices$resilience_capabilities,
                 matrix(1:6, 2, 3), matrix(2, 4, 4),
                 matrix(c(0, -1, 1, 0), 2, 2))) {
    reasons <- assumption_checks(A)$reason
    expect_false(any(grepl("NULL|NA_real_|stopifnot|\\bnrow\\b|list\\(", reasons)))
    expect_true(all(nchar(reasons) > 20))
  }
})

test_that("reasons read correctly when the count is one", {
  # These strings are read by users and the count is frequently 1 -- a single
  # blank row is the commonest way to fail a check. "1 factors are" is the
  # first thing a user would see on their own matrix.
  singular <- list(
    matrix(c(0, NA, 1, 0), 2, 2),
    matrix(c(0, -1, 1, 0), 2, 2),
    matrix(c(0, 3, 3, 1, 3,
             0, 0, 2, 0, 1,
             0, 0, 0, 2, 0,
             3, 1, 2, 0, 3,
             4, 1, 2, 1, 0), 5, 5, byrow = TRUE)
  )

  for (A in singular) {
    storage.mode(A) <- "double"
    reasons <- assumption_checks(A)$reason
    flagged <- reasons[grepl("^1 ", reasons)]
    expect_gt(length(flagged), 0)
    expect_false(any(grepl("^1 [a-z]+s ", flagged)))   # "1 factors", "1 entries"
    expect_false(any(grepl("^1 [a-z]+ are ", flagged)))
  }

  # And the plural still reads as a plural.
  plural <- assumption_checks(worked_matrices$resilience_capabilities)$reason
  expect_true(any(grepl("^4 factors are cut off", plural)))
})
