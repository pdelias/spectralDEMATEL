# The line that sets s, prominence and relation, and the leverage summary.
#
# What this file catches: the reduced form and sensitivity_matrix() disagreeing,
# a stabiliser found off the line that sets s, and a ranking of links returned
# without the condition number.
#
# The independent routes are the reduced form of the derivative, built from the
# eigenvectors of D and not of T, the paper's total_relation() run in a separate
# R process (data-raw/reference_values.R), and a two-factor system worked by
# hand.

fefo <- worked_matrices$fefo_stock_control

# The reduced form of the total derivative,
#   d lambda_max / d a_ij = (1 + lambda_max)^2 / s * (v_i u_j - mu_max * on_ij),
# with u the right and v the left dominant eigenvector scaled so that v'u = 1,
# and on_ij = 1 when link (i, j) lies in a row or a column whose total is s.
# Everything is computed from D = A / s. sensitivity_matrix() works from T.
reduced_form <- function(A) {
  s <- max(rowSums(A), colSums(A))
  D <- A / s
  e <- eigen(D)
  k <- which.max(Re(e$values))
  mu <- Re(e$values[k])
  u <- Re(e$vectors[, k])
  el <- eigen(t(D))
  v <- Re(el$vectors[, which.max(Re(el$values))])
  v <- v / sum(v * u)
  lam <- mu / (1 - mu)
  rs <- rowSums(A)
  cs <- colSums(A)
  on <- outer(abs(rs - s) < 1e-12, abs(cs - s) < 1e-12, "|")
  total <- (1 + lam)^2 / s * (outer(v, u) - mu * on)
  diag(total) <- NA
  weight <- outer(v, u)
  diag(weight) <- NA
  list(total = total, on = on, mu = mu, weight = weight)
}

leverage_fixtures <- function() {
  set.seed(19)
  out <- list(fefo = fefo, cycle3 = cycle3, defective = defective,
              two_by_two = two_by_two)
  for (n in c(4, 6, 9, 14)) out[[paste0("random", n)]] <- random_admissible(n)
  lapply(out, function(A) { attributes(A) <- list(dim = dim(A)); A })
}

# ---- normalising_line -------------------------------------------------------

test_that("normalising_line names the row or the column that sets s", {
  by_row <- matrix(c(0, 3, 3,
                     1, 0, 1,
                     1, 1, 0), 3, 3, byrow = TRUE)
  l <- normalising_line(by_row)
  expect_named(l, c("s", "rows", "cols", "tied"))
  expect_equal(l$s, 6)
  expect_identical(l$rows, 1L)
  expect_identical(l$cols, integer(0))
  expect_false(l$tied)

  l2 <- normalising_line(t(by_row))
  expect_equal(l2$s, 6)
  expect_identical(l2$rows, integer(0))
  expect_identical(l2$cols, 1L)
  expect_false(l2$tied)
})

test_that("normalising_line reports a tie between lines", {
  # Row 1 and column 3 both total 5.
  tie <- matrix(c(0, 2, 3,
                  1, 0, 2,
                  1, 1, 0), 3, 3, byrow = TRUE)
  l <- normalising_line(tie)
  expect_equal(l$s, 5)
  expect_identical(l$rows, 1L)
  expect_identical(l$cols, 3L)
  expect_true(l$tied)

  # Two rows total 3, and every column totals 2.
  two_rows <- matrix(c(0, 1, 1, 1,
                       1, 0, 1, 1,
                       1, 0, 0, 0,
                       0, 1, 0, 0), 4, 4, byrow = TRUE)
  l2 <- normalising_line(two_rows)
  expect_equal(l2$s, 3)
  expect_identical(l2$rows, c(1L, 2L))
  expect_identical(l2$cols, integer(0))
  expect_true(l2$tied)
})

test_that("normalising_line agrees with dematel() on s", {
  for (A in leverage_fixtures()) {
    expect_equal(normalising_line(A)$s, dematel(A)$s)
  }
})

test_that("uniform totals still have a line, and every line ties", {
  # T is undefined for this matrix, yet which lines set s is still a fact about
  # it, and the panel screening needs it.
  A <- matrix(2, 4, 4)
  diag(A) <- 0
  expect_null(dematel(A))
  l <- normalising_line(A)
  expect_equal(l$s, 6)
  expect_identical(l$rows, 1:4)
  expect_identical(l$cols, 1:4)
  expect_true(l$tied)
})

# ---- sensitivity_matrix and the reduced form ---------------------------------

test_that("sensitivity_matrix keeps its three components and gains the line", {
  s <- sensitivity_matrix(fefo)
  expect_named(s, c("total", "local", "normalization", "line"))
  expect_identical(s$line, normalising_line(fefo))
  expect_equal(s$total, s$local + s$normalization)
})

test_that("the total derivative equals the reduced form", {
  # Tolerance. The two routes use different eigendecompositions, of T and of
  # D. The derivative is of order (1 + lambda_max)^2 / s, so it is compared
  # relative to its largest entry. The defective fixture is included: its
  # repeated eigenvalue is the subdominant one, and the dominant pair that the
  # derivative uses is simple.
  for (nm in names(leverage_fixtures())) {
    A <- leverage_fixtures()[[nm]]
    s <- sensitivity_matrix(A)
    r <- reduced_form(A)
    scale <- max(abs(r$total), na.rm = TRUE)
    expect_lt(max(abs(s$total - r$total), na.rm = TRUE) / scale, 1e-9)
  }
})

test_that("no stabiliser lies off the line that normalising_line reports", {
  seen <- 0
  for (nm in names(leverage_fixtures())) {
    A <- leverage_fixtures()[[nm]]
    s <- sensitivity_matrix(A)
    n <- nrow(A)
    on <- matrix(FALSE, n, n)
    on[s$line$rows, ] <- TRUE
    on[, s$line$cols] <- TRUE
    neg <- which(s$total < 0, arr.ind = TRUE)
    seen <- seen + nrow(neg)
    expect_true(all(on[neg]), info = nm)
    # Off the line the normalisation effect is exactly zero.
    off_line <- !on & row(A) != col(A)
    expect_true(all(s$normalization[off_line] == 0), info = nm)
  }
  # The fixtures do hold stabilisers, so the assertion above is not vacuous.
  expect_gt(seen, 10)
})

test_that("a link is a stabiliser exactly when it is on the line and v_i u_j < mu_max", {
  for (nm in names(leverage_fixtures())) {
    A <- leverage_fixtures()[[nm]]
    s <- sensitivity_matrix(A)
    r <- reduced_form(A)
    closed <- r$on & r$weight < r$mu
    diag(closed) <- NA
    # Leave out links within rounding of the boundary.
    clear <- abs(r$weight - r$mu) > 1e-9 & !is.na(closed)
    expect_identical((s$total < 0)[clear], closed[clear], info = nm)
  }
})

# ---- prominence_relation ------------------------------------------------------

test_that("prominence and relation match a system worked by hand", {
  # A = [0 2; 1 0], s = 2, D = [0 1; 1/2 0].
  # (I - D)^-1 = [2 2; 1 2], so T = D (I - D)^-1 = [1 2; 1 1].
  # Row totals 3 and 2, column totals 2 and 3.
  p <- prominence_relation(two_by_two)
  expect_s3_class(p, "data.frame")
  expect_named(p, c("factor", "dispatch", "receipt", "prominence", "relation",
                    "role"))
  expect_equal(p$dispatch,   c(3, 2))
  expect_equal(p$receipt,    c(2, 3))
  expect_equal(p$prominence, c(5, 5))
  expect_equal(p$relation,   c(1, -1))
  expect_identical(p$role,   c("net cause", "net effect"))
})

test_that("prominence and relation match the paper's totals on the FEFO matrix", {
  # rowSums and colSums of total_relation(dematel_normalise(A, "keep")) in the
  # paper's code.
  dispatch <- c(4.7022389739652048, 4.6296338926487053, 5.35890428602473,
                5.555185736108788, 5.7667425219803414, 3.7835940767330882,
                4.8186770385019697, 3.3657579323994331, 4.9494978511905803,
                4.473206238336739, 4.6311621587596106, 5.0961131181785886,
                5.1386693132520236)
  receipt <- c(5.6312682246054555, 4.6896878103624289, 4.5930830748271418,
               4.2267242427029972, 5.1639197259097935, 3.4361849760151477,
               4.8915579156774394, 5.1387131421843488, 3.8534338349611446,
               4.1403048158041731, 5.2830956070592459, 5.5365009224579138,
               5.6849088455125729)
  p <- prominence_relation(fefo)
  expect_equal(p$dispatch, dispatch, tolerance = 1e-10)
  expect_equal(p$receipt,  receipt,  tolerance = 1e-10)
  expect_equal(p$prominence, dispatch + receipt, tolerance = 1e-10)
  expect_equal(p$relation,   dispatch - receipt, tolerance = 1e-10)
  expect_identical(p$role, ifelse(dispatch > receipt, "net cause", "net effect"))
})

test_that("relation sums to zero and the factors keep their names", {
  set.seed(2)
  A <- random_admissible(7)
  dimnames(A) <- list(paste0("F", 1:7), paste0("F", 1:7))
  p <- prominence_relation(A)
  expect_equal(sum(p$relation), 0, tolerance = 1e-10)
  expect_equal(sum(p$dispatch), sum(p$receipt))
  expect_identical(p$factor, paste0("F", 1:7))

  # Without names the factors are numbered.
  expect_identical(prominence_relation(two_by_two)$factor, c("1", "2"))
})

test_that("symmetric ratings give a relation of zero up to rounding", {
  # Symmetric ratings give equal row and column totals of T. The role of such
  # a factor is decided by rounding, which is why the documentation says to
  # read the sign with the size. The NA role of an exact zero is not asserted
  # here: no admissible fixture was found that produces one in floating point.
  A <- matrix(c(0, 1, 2,
                1, 0, 3,
                2, 3, 0), 3, 3, byrow = TRUE)
  p <- prominence_relation(A)
  expect_equal(p$relation, c(0, 0, 0), tolerance = 1e-12)
  expect_equal(p$dispatch, p$receipt, tolerance = 1e-12)
})

# ---- leverage_summary ---------------------------------------------------------

test_that("leverage_summary returns the ranking and the condition number together", {
  l <- leverage_summary(fefo)
  expect_named(l, c("line", "stabilisers", "amplifiers", "ev_condition",
                    "mu_max", "strongly_connected", "engine_version"))
  expect_true(l$strongly_connected)
  expect_identical(l$line, normalising_line(fefo))
  expect_equal(l$ev_condition, spectral_diagnostics(fefo)$ev_condition)
  expect_equal(l$mu_max, spectral_diagnostics(fefo)$mu_max)
  expect_identical(l$engine_version,
                   as.character(utils::packageVersion("spectralDEMATEL")))
  expect_named(l$stabilisers, c("from", "to", "rating", "total", "local",
                                "normalization"))
  expect_named(l$amplifiers, c("from", "to", "rating", "weight", "total"))
})

test_that("the stabilisers and amplifiers are those of the reduced form", {
  for (nm in names(leverage_fixtures())) {
    A <- leverage_fixtures()[[nm]]
    n <- nrow(A)
    l <- leverage_summary(A)
    r <- reduced_form(A)

    neg <- which(r$total < 0, arr.ind = TRUE)
    pos <- which(r$total > 0, arr.ind = TRUE)
    expect_equal(nrow(l$stabilisers), nrow(neg), info = nm)
    expect_equal(nrow(l$amplifiers),  nrow(pos), info = nm)
    expect_equal(nrow(l$stabilisers) + nrow(l$amplifiers), n * (n - 1), info = nm)

    key <- function(from, to) paste(from, to)
    expect_setequal(key(l$stabilisers$from, l$stabilisers$to),
                    key(neg[, 1], neg[, 2]))

    # Stabilisers: every one on the line, most negative first.
    if (nrow(l$stabilisers) > 0) {
      expect_true(all(l$stabilisers$from %in% l$line$rows |
                        l$stabilisers$to %in% l$line$cols), info = nm)
      expect_true(all(l$stabilisers$total < 0))
      expect_false(is.unsorted(l$stabilisers$total))
      expect_equal(l$stabilisers$total,
                   l$stabilisers$local + l$stabilisers$normalization)
    }

    # Amplifiers: ranked by v_i u_j with v'u = 1, largest first.
    a <- l$amplifiers
    expect_false(is.unsorted(rev(a$weight)))
    scale <- max(abs(r$weight), na.rm = TRUE)
    expect_lt(max(abs(a$weight - r$weight[cbind(a$from, a$to)])) / scale, 1e-9)
    expect_equal(a$rating, A[cbind(a$from, a$to)])
    expect_true(all(a$total > 0))
  }
})

test_that("the two-factor system has one stabiliser and one amplifier", {
  # A = [0 2; 1 0]. Row 1 and column 2 both total 2 and set s, so link (1, 2)
  # is on the line and link (2, 1) is off it. With D = [0 1; 1/2 0],
  # mu_max = sqrt(1/2). The right vector is (sqrt(2), 1) and the left
  # (1, sqrt(2)), up to scale, so v'u = 2 sqrt(2) and the weights are
  # v_1 u_2 = 1 / (2 sqrt(2)) = 0.354 and v_2 u_1 = 2 / (2 sqrt(2)) = 0.707.
  # Link (1, 2): on the line, weight 0.354 < mu_max = 0.707, a stabiliser.
  # Link (2, 1): off the line, an amplifier with weight 0.707.
  l <- leverage_summary(two_by_two)
  expect_true(l$line$tied)
  expect_equal(nrow(l$stabilisers), 1)
  expect_equal(c(l$stabilisers$from, l$stabilisers$to), c(1, 2))
  expect_equal(nrow(l$amplifiers), 1)
  expect_equal(c(l$amplifiers$from, l$amplifiers$to), c(2, 1))
  expect_equal(l$amplifiers$weight, 1 / sqrt(2), tolerance = 1e-12)
  expect_equal(l$mu_max, sqrt(1 / 2), tolerance = 1e-12)
})

test_that("a system that is not strongly connected gets no rows made of noise", {
  # Three matrices found by the review of 0.7.0. Each is accepted by dematel()
  # and is not strongly connected. Before the guard, the first returned tables
  # filled with NA, the second listed as stabilisers four links off the line
  # whose derivative was a rounding error of 1e-17, and the third ranked links
  # of a repeated dominant eigenvalue.
  cases <- list(
    one_link = matrix(c(0, 3, 0, 0), 2, 2, byrow = TRUE),
    zeros_in_profile = rbind(c(0, 0, 0, 3), c(3, 3, 3, 1), c(3, 2, 1, 1),
                             c(1, 0, 0, 1)),
    two_blocks = block_diagonal * c(1, 2))

  for (nm in names(cases)) {
    A <- cases[[nm]]
    expect_false(is_irreducible(A), info = nm)
    expect_silent(l <- leverage_summary(A))
    expect_false(l$strongly_connected, info = nm)
    for (tab in list(l$stabilisers, l$amplifiers)) {
      expect_false(anyNA(tab), info = nm)
    }
    expect_true(all(l$stabilisers$from %in% l$line$rows |
                      l$stabilisers$to %in% l$line$cols), info = nm)
    expect_true(all(l$stabilisers$total < 0), info = nm)
    expect_true(all(l$amplifiers$total > 0), info = nm)
  }
  # The second case: one row sets s, and no stabiliser lies outside it.
  l <- leverage_summary(cases$zeros_in_profile)
  expect_identical(l$line$rows, 2L)
  expect_true(all(l$stabilisers$from == 2L))
})

test_that("normalising_line applies an absolute tolerance, as documented", {
  # Row 1 totals 0.1 + 0.2, which is 0.30000000000000004 in floating point, and
  # column 1 totals 0.3. The rows are considered because their maximum is the
  # larger by exact comparison, so the column is not reported although it ties
  # in exact arithmetic. This is the rule of 0.6.0 and it is pinned here so
  # that it does not change unnoticed.
  A <- rbind(c(0, 0.1, 0.2), c(0.3, 0, 0), c(0, 0, 0))
  l <- normalising_line(A)
  expect_identical(l$rows, c(1L, 2L))
  expect_identical(l$cols, integer(0))
})
