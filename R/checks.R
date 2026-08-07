# The assumption checks, returned as data rather than raised as conditions.
#
# The paper's Section 2 conditions become values the engine returns. Three
# consumers need this shape and warnings serve none of them: an interface
# renders the checks as a panel, a batch job tabulates them across a hundred
# matrices, and a diagnosis card prints them as provenance. A warning() can only
# be caught and re-parsed, which is how message text becomes an accidental API.

# Thresholds for the two continuous checks. Neither is fitted; both are
# recommendations, chosen against measurements rather than invented, and both
# are documented in assumption_checks() so an interface can say so.
#
#   coupling_warn      The corpus of 117 published systems has a 99th percentile
#                      coupling of 0.961 and a maximum of 0.965, so 0.95 flags
#                      roughly the top 3% -- systems close enough to criticality
#                      that the multiplier is very sensitive to small changes.
#
#   ev_condition_warn  Measured across 300 random admissible matrices of 5 to 30
#                      factors, ev_condition ran from 1.0 to 2.0. A symmetric
#                      hub gives exactly 1, an asymmetric hub 3.2, and a long
#                      directed chain 28 to 33. Chains and deep hierarchies are
#                      where a first-order derivative stops being informative,
#                      and 5 sits in the empty gap between the two regimes.
CHECK_THRESHOLDS <- list(
  coupling_warn     = 0.95,
  ev_condition_warn = 5
)

# The checks in the order they are reported, which is also the order in which
# assumption_checks() emits them and the order skip_remaining() indexes into.
# Fixed, so a batch job can rbind a hundred of these and get a rectangular
# table. Change one and you must change the other.
CHECK_IDS <- c(
  "input_is_matrix",          # 1  structural: nothing else runs without these
  "square",                   # 2
  "finite",                   # 3
  "nonnegative",              # 4
  "has_influence",            # 5
  "zero_diagonal",            # 6  convention
  "strong_connectivity",      # 7  the graph
  "connectivity_margin",      # 8
  "totals_vary",              # 9  the normalisation
  "coupling_margin",          # 10 the spectrum
  "sensitivity_conditioning"  # 11
)

#' Check a matrix against the assumptions the diagnostics rest on
#'
#' Every condition is returned as a row of data: a verdict, an identifier stable
#' enough to branch on, a reason written for a user rather than a developer, and
#' where possible the factors at fault. Nothing here stops, warns or prints.
#'
#' @param A The matrix to check. Anything may be passed; a non-matrix fails the
#'   first check rather than raising an error.
#' @param type Whether `A` is a direct influence matrix (`"A"`) or an
#'   already-published total-relation matrix (`"T"`). Two checks do not apply to
#'   a total-relation matrix and are reported as `skipped`.
#'
#' @return A data frame of class `spectral_checks`, one row per check, in a
#'   fixed order, with columns:
#'   \describe{
#'     \item{`check`}{Stable identifier. Branch on this, never on `reason`.}
#'     \item{`verdict`}{One of `"pass"`, `"warn"`, `"fail"`, `"skipped"`.}
#'     \item{`reason`}{Plain-language explanation, and what to do about it.}
#'     \item{`value`}{The number behind the verdict, or `NA` where there is
#'       none.}
#'     \item{`factors`}{List column of integer indices of the factors at fault,
#'       empty where the check is not about particular factors.}
#'   }
#'
#' @section The four verdicts:
#' `pass` and `fail` are the assumption holding or not. `warn` is the case a
#' binary verdict cannot express: the matrix is in scope, but it is close enough
#' to a boundary that the reading is worth less than it looks. Strong
#' connectivity is yes or no, but a graph one edge from disconnection behaves
#' differently from a dense one, and numerical conditioning degrades long before
#' any check trips.
#'
#' `skipped` means the check could not be evaluated because something it depends
#' on failed first — there is no point testing connectivity on a matrix that is
#' not square. It is not a pass, and an interface must not render it as one.
#'
#' @section A failure is not a refusal:
#' The engine computes diagnostics for a matrix that fails these checks, and
#' [spectral_diagnostics()] returns the checks alongside the numbers so the two
#' cannot be separated. One of the two matrices in [worked_matrices] fails
#' strong connectivity and its diagnostics are published. Deciding what to do
#' about a failure belongs to whoever is calling this.
#'
#' @examples
#' checks <- assumption_checks(worked_matrices$resilience_capabilities)
#' checks[checks$verdict != "pass", c("check", "verdict", "reason")]
#'
#' # Which factor is at fault:
#' checks$factors[checks$check == "strong_connectivity"][[1]]
#'
#' @seealso [spectral_diagnostics()], which returns these alongside the numbers.
#' @export
assumption_checks <- function(A, type = c("A", "T")) {
  type <- match.arg(type)
  rows <- list()
  add <- function(...) rows[[length(rows) + 1L]] <<- check_row(...)

  # ---- structural: everything else depends on these --------------------------

  ok_matrix <- is.matrix(A) && is.numeric(A)
  add("input_is_matrix", if (ok_matrix) "pass" else "fail",
      if (ok_matrix) "The input is a numeric matrix."
      else paste("The input is not a numeric matrix. A data frame read from a",
                 "spreadsheet usually needs as.matrix(), and a leading column",
                 "of factor names must be removed first."))

  if (!ok_matrix) return(skip_remaining(rows, from = 2L))

  ok_square <- nrow(A) == ncol(A)
  add("square", if (ok_square) "pass" else "fail",
      if (ok_square) sprintf("The matrix is square, with %d factors.", nrow(A))
      else sprintf(paste("The matrix has %d rows and %d columns. An influence",
                         "matrix must be square, one row and one column per",
                         "factor; a header row or an index column is the usual",
                         "cause."), nrow(A), ncol(A)),
      value = nrow(A))

  if (!ok_square) return(skip_remaining(rows, from = 3L))

  bad_finite <- which(!is.finite(A), arr.ind = TRUE)
  ok_finite <- nrow(bad_finite) == 0L
  add("finite", if (ok_finite) "pass" else "fail",
      if (ok_finite) "Every entry is a finite number."
      else sprintf(paste("%s missing or infinite. Blank cells in a",
                         "spreadsheet read as missing; an unrated pair should",
                         "be entered as 0."),
                   count_phrase(nrow(bad_finite), "entry", "is", "are")),
      value = nrow(bad_finite),
      factors = sort(unique(as.integer(bad_finite[, "row"]))))

  bad_neg <- if (ok_finite) which(A < 0, arr.ind = TRUE) else
    matrix(integer(0), 0, 2, dimnames = list(NULL, c("row", "col")))
  ok_neg <- ok_finite && nrow(bad_neg) == 0L
  add("nonnegative", if (!ok_finite) "skipped" else if (ok_neg) "pass" else "fail",
      if (!ok_finite) "Not evaluated: the matrix contains missing or infinite entries."
      else if (ok_neg) "Every entry is zero or positive."
      else sprintf(paste("%s negative. DEMATEL ratings measure the",
                         "strength of an influence, not its direction, so a",
                         "negative entry usually means a sign convention from",
                         "another method has been carried over."),
                   count_phrase(nrow(bad_neg), "entry", "is", "are")),
      value = nrow(bad_neg),
      factors = sort(unique(as.integer(bad_neg[, "row"]))))

  structural_ok <- ok_finite && ok_neg
  if (!structural_ok) return(skip_remaining(rows, from = 5L))

  ok_any <- any(A > 0)
  add("has_influence", if (ok_any) "pass" else "fail",
      if (ok_any) "The matrix records at least one influence."
      else "Every entry is zero, so there is no system to diagnose.",
      value = sum(A > 0))

  if (!ok_any) return(skip_remaining(rows, from = 6L))

  # ---- conventions and the graph --------------------------------------------

  diag_nonzero <- which(diag(A) != 0)
  add("zero_diagonal", if (length(diag_nonzero) == 0L) "pass" else "warn",
      if (length(diag_nonzero) == 0L)
        "The diagonal is zero, as the convention expects."
      else sprintf(paste("%s a non-zero self-influence. This is",
                         "admissible and the diagnostics are computed as usual,",
                         "but self-influence inflates coupling, so the value is",
                         "not comparable with studies that set the diagonal to",
                         "zero."), count_phrase(length(diag_nonzero), "factor", "carries", "carry")),
      value = length(diag_nonzero),
      factors = as.integer(diag_nonzero))

  P <- reachability(A)
  connected <- all(P > 0)
  stranded <- stranded_factors(P)
  add("strong_connectivity", if (connected) "pass" else "fail",
      if (connected)
        "Every factor can reach every other factor, so assumption A2 holds."
      else sprintf(paste("%s cut off from the main body of the",
                         "system, so assumption A2 fails and the entry profile",
                         "is no longer unique. This is usually a coding slip in",
                         "one row rather than a property of the system: check",
                         "whether those factors were left blank."),
                   count_phrase(length(stranded), "factor", "is", "are")),
      value = length(stranded),
      factors = as.integer(stranded))

  M <- (A > 0) * 1
  fragile <- sort(unique(c(which(rowSums(M) == 1), which(colSums(M) == 1))))
  add("connectivity_margin",
      if (!connected) "skipped" else if (length(fragile) == 0L) "pass" else "warn",
      if (!connected)
        "Not evaluated: the influence graph is already disconnected."
      else if (length(fragile) == 0L)
        "No factor depends on a single link, so connectivity is not fragile."
      else sprintf(paste("%s on a single incoming or outgoing",
                         "link. Strong connectivity holds, but removing one",
                         "rating would break it, so the diagnosis is less",
                         "robust than a denser matrix of the same size."),
                   count_phrase(length(fragile), "factor", "hangs", "hang")),
      # Passed unconditionally: check_row() drops both when the verdict is
      # skipped, which is where that invariant is enforced.
      value = length(fragile),
      factors = as.integer(fragile))

  # ---- the normalisation and the spectrum -----------------------------------

  if (type == "T") {
    add("totals_vary", "skipped",
        paste("Not applicable: a total-relation matrix is not normalised, so",
              "the uniform-totals pathology cannot arise."))
    m <- NULL
    d <- try_diagnostics(A, "T")
  } else {
    m <- dematel(A)
    add("totals_vary", if (!is.null(m)) "pass" else "fail",
        if (!is.null(m))
          "Row and column totals vary across factors, so the total-relation matrix exists."
        else paste("Every row and column total is the same. The normalised",
                   "matrix then has spectral radius exactly 1, the total-relation",
                   "matrix is undefined, and no diagnostic can be computed. This",
                   "is a property of the input, not a mistake: a matrix where",
                   "every pair was rated identically will always meet it."))
    d <- if (is.null(m)) NULL else try_diagnostics(A, "A")
  }

  if (is.null(d)) return(skip_remaining(rows, from = length(rows) + 1L))

  near_critical <- d$mu_max >= CHECK_THRESHOLDS$coupling_warn
  add("coupling_margin", if (near_critical) "warn" else "pass",
      if (!near_critical)
        sprintf(paste("Coupling is %.3f, comfortably below criticality; the",
                      "total-effect multiplier of %.2f is stable under small",
                      "changes to the ratings."), d$mu_max, d$multiplier)
      else sprintf(paste("Coupling is %.3f, within %.3f of criticality. The",
                         "multiplier of %.2f grows very fast in this range, so a",
                         "small change to one rating can move it a long way.",
                         "Treat the magnitude as indicative rather than exact."),
                   d$mu_max, 1 - d$mu_max, d$multiplier),
      value = d$mu_max)

  ill_conditioned <- d$ev_condition >= CHECK_THRESHOLDS$ev_condition_warn
  add("sensitivity_conditioning", if (ill_conditioned) "warn" else "pass",
      if (!ill_conditioned)
        sprintf(paste("The eigenvalue condition number is %.2f, so the per-link",
                      "sensitivity estimates can be read at face value."),
                d$ev_condition)
      else sprintf(paste("The eigenvalue condition number is %.1f. First-order",
                         "sensitivity estimates lose roughly that factor of",
                         "accuracy, so the per-link ranking is indicative at",
                         "best. Long chains and deep hierarchies produce this."),
                   d$ev_condition),
      value = d$ev_condition)

  checks_to_df(rows)
}

# ---- internals -------------------------------------------------------------

#' "1 factor is" rather than "1 factors are".
#'
#' These strings are read by users, and the counts they interpolate are
#' frequently 1 -- a single blank row is the commonest way to fail a check.
#' @noRd
count_phrase <- function(n, noun, verb_singular, verb_plural) {
  sprintf("%d %s%s %s", n, noun, if (n == 1) "" else "s",
          if (n == 1) verb_singular else verb_plural)
}

#' One row of the checks table.
#'
#' Enforces the invariant that a `skipped` check carries neither a value nor
#' factors. Two call sites got this wrong independently -- one reported the
#' factors hanging on a single link when the graph was already disconnected,
#' the other reported a count of zero negative entries on a matrix whose
#' entries had never been examined. A number or a list of names beside a
#' verdict that was never reached reads as evidence for it. Enforcing it here
#' means no future call site can reintroduce the problem.
#' @noRd
check_row <- function(check, verdict, reason, value = NA_real_,
                      factors = integer(0)) {
  if (identical(verdict, "skipped")) {
    value <- NA_real_
    factors <- integer(0)
  }
  list(check = check, verdict = verdict, reason = reason,
       value = as.numeric(value), factors = as.integer(factors))
}

#' Assemble rows into the returned data frame, keeping `factors` as a list
#' column so a caller gets indices rather than a string to re-parse.
#' @noRd
checks_to_df <- function(rows) {
  out <- data.frame(
    check   = vapply(rows, `[[`, character(1), "check"),
    verdict = vapply(rows, `[[`, character(1), "verdict"),
    reason  = vapply(rows, `[[`, character(1), "reason"),
    value   = vapply(rows, `[[`, numeric(1),   "value"),
    stringsAsFactors = FALSE
  )
  out$factors <- lapply(rows, `[[`, "factors")
  rownames(out) <- NULL
  class(out) <- c("spectral_checks", "data.frame")
  out
}

#' Fill the remaining checks with `skipped` so the table is always rectangular
#' and a batch job can rbind a hundred of them.
#' @noRd
skip_remaining <- function(rows, from) {
  for (id in CHECK_IDS[seq(from, length(CHECK_IDS))]) {
    rows[[length(rows) + 1L]] <- check_row(
      id, "skipped",
      "Not evaluated: an earlier check the matrix had to pass first did not.")
  }
  checks_to_df(rows)
}

#' Boolean reachability closure of the influence graph.
#'
#' The single implementation behind both [is_irreducible()] and the
#' `strong_connectivity` check, so the verdict and the named factors can never
#' disagree with each other.
#' @noRd
reachability <- function(A) {
  n <- nrow(A)
  M <- (A > 0) * 1
  R <- diag(n) + M
  P <- R
  for (k in seq_len(n)) P <- (P %*% R > 0) * 1
  P
}

#' Factors cut off from the largest mutually-reachable group.
#'
#' The obvious definition -- every factor that cannot reach, or be reached by,
#' all the others -- names the whole matrix as soon as one factor is stranded,
#' because nobody can reach that one either. Useless to a user looking for the
#' row they mistyped. This instead partitions into strongly connected components
#' (i and j share one exactly when each can reach the other) and reports what
#' falls outside the largest. On the resilience matrix that is factor 14 alone,
#' which is the row that is actually blank.
#'
#' Returns `integer(0)` for a strongly connected graph. Where two components are
#' the same size the choice of "largest" is arbitrary, and the reason text says
#' "main body" rather than claiming more than that.
#' @noRd
stranded_factors <- function(P) {
  mutual <- (P > 0) & (t(P) > 0)
  main <- which.max(rowSums(mutual))
  which(!mutual[main, ])
}

#' Diagnostics without the checks, returning NULL rather than propagating a
#' numerical failure out of a function whose whole job is to report problems.
#' @noRd
try_diagnostics <- function(A, type) {
  out <- tryCatch(spectral_diagnostics(A, type = type, checks = FALSE),
                  error = function(e) NULL)
  if (is.null(out) || !is.finite(out$mu_max)) NULL else out
}
