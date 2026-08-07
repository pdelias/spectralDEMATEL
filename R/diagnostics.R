# The reported diagnostic set for one system.
#
# The arithmetic in diagnostics_core() is the paper's R/spectral_dematel.R
# unchanged, and that file remains the single source of truth for every
# definition. Do not restate a formula here from memory or from an article's
# prose.

#' Spectral diagnostics for one DEMATEL system
#'
#' Computes the full diagnostic set from a single influence matrix: how strongly
#' the system amplifies its own influence, whether one propagation pattern
#' dominates, how concentrated the entry points are, where influence enters and
#' where it accumulates, and how far a first-order sensitivity estimate can be
#' trusted.
#'
#' @param A A square numeric matrix. Either a direct influence matrix
#'   (`type = "A"`, the default) or an already-published total-relation matrix
#'   (`type = "T"`).
#' @param type Which matrix `A` is. `"A"` normalises and builds `T` through
#'   [dematel()]. `"T"` takes `A` as the total-relation matrix and recovers the
#'   spectrum of `D` through the inverse Moebius map \eqn{\mu = \lambda/(1 +
#'   \lambda)}.
#' @param checks Whether to run [assumption_checks()] and attach the result.
#'   Defaults to `TRUE` and should stay that way for anything a person will
#'   read. `FALSE` returns the numbers alone, or `NULL` if they cannot be
#'   computed, and exists for inner loops such as [surrogate_ensemble()].
#'
#' @return With `checks = TRUE`, a named list that always has the same shape,
#'   whatever was passed in. With `checks = FALSE`, the same list without
#'   `checks` and `engine_version`, or `NULL` when the matrix is inadmissible.
#'
#'   \describe{
#'     \item{`n`}{Number of factors.}
#'     \item{`mu_max`}{**Coupling**, in \eqn{(0, 1)}. The spectral radius of the
#'       normalised matrix, so the distance to criticality. Higher means more
#'       strongly coupled.}
#'     \item{`multiplier`}{Total-to-direct multiplier, \eqn{1/(1 - \mu_{max})}.
#'       Always greater than 1.}
#'     \item{`lambda_max`}{Dominant eigenvalue of the total-relation matrix.}
#'     \item{`indirect_dominant`}{`TRUE` when indirect effects exceed direct
#'       ones, which happens exactly when `lambda_max > 1`, equivalently when
#'       `mu_max > 0.5`.}
#'     \item{`dominance`}{**Mode dominance**, \eqn{|\lambda_2|/\lambda_{max}} in
#'       \eqn{[0, 1)}. Low means a single propagation mode; high means a second
#'       mode competes with the first.}
#'     \item{`hierarchy_sd`}{Hierarchy as the standard deviation of the
#'       unit-norm right eigenvector. **HIGH means concentrated.**}
#'     \item{`hierarchy_pr`}{Hierarchy as the participation ratio divided by
#'       `n`, in \eqn{(0, 1]}. **LOW means concentrated** — the opposite
#'       direction to the other two.}
#'     \item{`hierarchy_gini`}{Hierarchy as the Gini coefficient of the same
#'       vector. **HIGH means concentrated.** Size-free.}
#'     \item{`entry_points`}{Right eigenvector, unit 2-norm. Weights each factor
#'       by the influence it injects: where to apply pressure.}
#'     \item{`accumulation`}{Left eigenvector, unit 2-norm. Weights each factor
#'       by what accumulates there: where effects land.}
#'     \item{`ev_condition`}{Eigenvalue condition number, \eqn{1/|v'u|}. At least
#'       1 by construction. Bounds how far the first-order sensitivity estimates
#'       of [sensitivity_matrix()] can be trusted; large values mean the
#'       derivative is locally uninformative.}
#'     \item{`checks`}{The [assumption_checks()] table for this matrix.}
#'     \item{`engine_version`}{The package version that produced the diagnosis.}
#'   }
#'
#' @section The numbers never travel without their checks:
#' A diagnosis computed on a matrix that failed a check has to carry that fact
#' with it, so `checks` is part of the returned object rather than something the
#' caller was told separately and may not have kept. When the matrix is
#' inadmissible every numeric field is `NA` and `checks` says why — the shape
#' does not change, so a batch job over a hundred matrices gets a rectangular
#' result and never a failure.
#'
#' `engine_version` is recorded because a metric definition has moved here
#' before. When one moves again, old diagnoses must remain interpretable rather
#' than silently mean something new.
#'
#' @section The three hierarchy readings run in different directions:
#' `hierarchy_sd` and `hierarchy_gini` are high when influence enters at a few
#' factors. `hierarchy_pr` is **low** when influence enters at a few factors.
#' This project has inverted the direction once already, in code that had been
#' reviewed, and the mistake produced entirely plausible output. Anything
#' displaying one of these must state which reading it used and which way it
#' runs.
#'
#' `hierarchy_sd` is the reading the source paper's corpus analysis uses
#' throughout, so it is the one behind every corpus-level result including the
#' coupling-hierarchy trade-off.
#'
#' @section Mode dominance uses the largest modulus:
#' \eqn{\lambda_2} is the largest-modulus eigenvalue below the dominant one, not
#' the second largest real part. For a non-negative `T` the negative end of the
#' spectrum usually carries the larger modulus, so taking the real part
#' understates the ratio nearly everywhere — while passing every plausibility
#' check.
#'
#' @seealso [assumption_checks()] for the checks alone,
#'   [sensitivity_matrix()] for where to intervene,
#'   [surrogate_ensemble()] for a baseline to read these numbers against.
#'
#' @examples
#' A <- matrix(c(0, 2, 0,
#'               0, 0, 3,
#'               4, 0, 0), 3, 3, byrow = TRUE)
#' d <- spectral_diagnostics(A)
#' d$mu_max
#' d$multiplier
#' d$indirect_dominant
#'
#' # Anything that failed, and why:
#' subset(d$checks, verdict != "pass", c(check, verdict))
#'
#' # Inadmissible input keeps the shape and explains itself:
#' bad <- spectral_diagnostics(matrix(2, 4, 4))
#' bad$mu_max
#' bad$checks$reason[bad$checks$check == "totals_vary"]
#'
#' @export
spectral_diagnostics <- function(A, type = c("A", "T"), checks = TRUE) {
  type <- match.arg(type)
  core <- diagnostics_core(A, type)

  if (!checks) return(core)

  ck <- assumption_checks(A, type = type)
  if (is.null(core)) core <- na_diagnostics(A)
  c(core, list(checks = ck, engine_version = engine_version()))
}

#' The arithmetic, ported unchanged from the paper's implementation.
#' Returns NULL for any matrix it cannot process.
#' @noRd
diagnostics_core <- function(A, type) {
  if (type == "T") {
    if (!is.matrix(A) || !is.numeric(A) || nrow(A) != ncol(A)) return(NULL)
    if (!all(is.finite(A))) return(NULL)
    # Recover the spectrum of D from a published total relation matrix.
    eT <- eigen(A)
    lam <- eT$values
    k <- which.max(Re(lam))
    lambda_max <- Re(lam[k])
    mu <- lam / (1 + lam)
    mu_max <- Re(mu[k])
    vecs <- eT$vectors
    lam2 <- sort(Mod(lam[-k]), decreasing = TRUE)[1]
    u <- abs(Re(vecs[, k]))
    vT <- eigen(t(A))
    v <- abs(Re(vT$vectors[, which.max(Re(vT$values))]))
  } else {
    m <- dematel(A)
    if (is.null(m)) return(NULL)
    eT <- eigen(m$T)
    lam <- eT$values
    k <- which.max(Re(lam))
    lambda_max <- Re(lam[k])
    mu_max <- lambda_max / (1 + lambda_max)
    lam2 <- sort(Mod(lam[-k]), decreasing = TRUE)[1]
    u <- abs(Re(eT$vectors[, k]))
    vT <- eigen(t(m$T))
    v <- abs(Re(vT$vectors[, which.max(Re(vT$values))]))
  }
  if (sum(u^2) == 0 || sum(v^2) == 0) return(NULL)
  u <- u / sqrt(sum(u^2))
  v <- v / sqrt(sum(v^2))
  n <- length(u)

  list(
    n                 = n,
    mu_max            = mu_max,                       # coupling: distance to criticality, in (0,1)
    multiplier        = 1 / (1 - mu_max),             # total-to-direct multiplier (= 1 + lambda_max)
    lambda_max        = lambda_max,
    indirect_dominant = lambda_max > 1,               # indirect effects exceed direct effects
    dominance         = lam2 / lambda_max,            # |lambda_2|/lambda_max, in [0,1); low = single mode
    # Three readings of the same idea, in DIFFERENT directions. hierarchy_sd is
    # the one the source corpus pipeline reports.
    hierarchy_sd      = stats::sd(u),                 # SD of u (n-1 denom); HIGH = concentrated
    hierarchy_pr      = (sum(u^2)^2 / sum(u^4)) / n,  # participation ratio / n, in (0,1]; LOW = concentrated
    hierarchy_gini    = gini(u),                      # size-free; HIGH = concentrated
    entry_points      = u,                            # right eigenvector: where to inject pressure
    accumulation      = v,                            # left eigenvector: where effects land
    ev_condition      = 1 / abs(sum(v * u))           # reliability of sensitivity estimates (>= 1)
  )
}

#' The same shape, with every number missing, for a matrix whose diagnostics
#' cannot be computed. Keeps a batch result rectangular.
#' @noRd
na_diagnostics <- function(A) {
  n <- if (is.matrix(A) && is.numeric(A) && nrow(A) == ncol(A)) nrow(A) else NA_integer_
  vec <- if (is.na(n)) NA_real_ else rep(NA_real_, n)
  list(
    n                 = n,
    mu_max            = NA_real_,
    multiplier        = NA_real_,
    lambda_max        = NA_real_,
    indirect_dominant = NA,
    dominance         = NA_real_,
    hierarchy_sd      = NA_real_,
    hierarchy_pr      = NA_real_,
    hierarchy_gini    = NA_real_,
    entry_points      = vec,
    accumulation      = vec,
    ev_condition      = NA_real_
  )
}

#' The version stamped onto every diagnosis.
#' @noRd
engine_version <- function() {
  as.character(utils::packageVersion("spectralDEMATEL"))
}

#' Gini coefficient of a non-negative vector
#'
#' The size-free concentration reading behind `hierarchy_gini`. Uses the
#' rank-weighted form, which is \eqn{O(n \log n)} rather than the \eqn{O(n^2)}
#' mean-absolute-difference form and agrees with it exactly.
#'
#' @param x A non-negative numeric vector.
#' @return A scalar in \eqn{[0, 1)}. High means concentrated.
#' @noRd
gini <- function(x) {
  x <- sort(x)
  n <- length(x)
  sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x))
}
